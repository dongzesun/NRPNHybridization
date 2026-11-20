"""Ordered chunks of radiation work over one read-only, process-shareable history."""

from collections import deque
from concurrent.futures import wait
from contextlib import contextmanager
from itertools import islice
from pathlib import Path
from tempfile import TemporaryDirectory
import pickle
import struct

import numpy as np


@contextmanager
def shared_context(context):
    """Snapshot once; large NumPy buffers are mapped read-only by local workers."""
    with TemporaryDirectory(prefix="eccpn-radiation-") as name:
        path = Path(name)
        with (path / "arrays").open("wb") as arrays, (path / "metadata").open("wb") as metadata:
            def write_buffer(buffer):
                view = buffer.raw()
                arrays.write(struct.pack("=Q", view.nbytes))
                arrays.write(view)
                arrays.write(b"\0" * (-view.nbytes % 8))

            pickle.dump(context, metadata, protocol=5, buffer_callback=write_buffer)
        yield path


def _run_chunk(path, function, job):
    sequence, items, maximum_chunk = job
    arrays = path / "arrays"
    mapped = np.memmap(arrays, dtype=np.uint8, mode="r") if arrays.stat().st_size else b""
    view = memoryview(mapped)

    def buffers():
        offset = 0
        while offset < len(view):
            length, = struct.unpack_from("=Q", view, offset)
            offset += 8
            yield view[offset:offset + length]
            offset += length + (-length % 8)

    with (path / "metadata").open("rb") as metadata:
        context = pickle.load(metadata, buffers=buffers())
    if len(items) <= maximum_chunk:
        return function(context, items), None
    # Amortize decoding a long history across a worker range without sending
    # that range's entire output through the executor queue at once.
    spool = path / f"result-{sequence}"
    with spool.open("wb") as output:
        for start in range(0, len(items), maximum_chunk):
            result = function(context, items[start:start + maximum_chunk])
            pickle.dump(result, output, protocol=5)
    return None, spool


def ordered_chunks(executor, context, function, items, maximum_chunk=32):
    """Run consecutive chunks with at most one pending result per worker.

    A chunk owns its mutable caches, including when executed by a thread.
    The caller keeps ownership of the executor. On failure all work submitted
    here has finished or been cancelled before its mapped files are removed.
    """
    if not len(items):
        return
    if executor is None or len(items) == 1:
        for start in range(0, len(items), maximum_chunk):
            yield function(context, items[start:start + maximum_chunk])
        return
    workers = max(1, getattr(executor, "_max_workers", 1))
    size = min((len(items) + workers - 1) // workers,
               max(maximum_chunk, (len(items) + 3 * workers - 1) // (3 * workers)))
    chunks = ((start, items[start:start + size], maximum_chunk) for start in range(0, len(items), size))
    with shared_context(context) as path:
        pending = deque()
        try:
            for chunk in islice(chunks, workers):
                pending.append(executor.submit(_run_chunk, path, function, chunk))
            while pending:
                result, spool = pending[0].result()
                pending.popleft()
                chunk = next(chunks, None)
                if chunk is not None:
                    pending.append(executor.submit(_run_chunk, path, function, chunk))
                if spool is None:
                    yield result
                else:
                    with spool.open("rb") as output:
                        while output.peek(1):
                            yield pickle.load(output)
                    spool.unlink()
        finally:
            for future in pending:
                future.cancel()
            wait(pending)
