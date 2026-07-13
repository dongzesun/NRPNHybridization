"""Optional provenance recording for BMS frame-fixing diagnostics."""

from __future__ import annotations

import hashlib
import json
import os
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import quaternion


def _hash_arrays(*arrays):
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(np.asarray(array))
        digest.update(str(value.dtype).encode("ascii"))
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.view(np.uint8))
    return digest.hexdigest()


def waveform_fingerprint(waveform):
    """Return a stable numerical fingerprint and compact waveform metadata."""
    return {
        "sha256": _hash_arrays(waveform.t, waveform.data),
        "n_times": int(len(waveform.t)),
        "n_modes": int(waveform.data.shape[1]),
        "time_range": [float(waveform.t[0]), float(waveform.t[-1])],
        "ell_min": int(waveform.ell_min),
        "ell_max": int(waveform.ell_max),
    }


def abd_fingerprint(abd):
    """Return a stable full-field fingerprint and compact ABD metadata."""
    digest = hashlib.sha256()
    for name in ("t", "sigma", "psi0", "psi1", "psi2", "psi3", "psi4"):
        value = getattr(abd, name)
        array = np.ascontiguousarray(np.asarray(value))
        digest.update(name.encode("ascii"))
        digest.update(str(array.dtype).encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.view(np.uint8))
    return {
        "sha256": digest.hexdigest(),
        "n_times": int(len(abd.t)),
        "time_range": [float(abd.t[0]), float(abd.t[-1])],
        "ell_max": int(abd.ell_max),
    }


def transform_summary(transformation):
    """Serialize the replay-relevant parts of a scri BMS transformation."""
    return {
        "ell_max": int(transformation.ell_max),
        "order": list(transformation.order),
        "boost_velocity": np.asarray(transformation.boost_velocity, dtype=float).tolist(),
        "frame_rotation_wxyz": quaternion.as_float_array(
            transformation.frame_rotation
        ).tolist(),
        "supertranslation_real": np.asarray(
            transformation.supertranslation
        ).real.tolist(),
        "supertranslation_imag": np.asarray(
            transformation.supertranslation
        ).imag.tolist(),
    }


def optimizer_summary(result):
    """Serialize a scipy OptimizeResult without large Jacobian arrays."""
    return {
        "x": np.asarray(result.x, dtype=float).tolist(),
        "cost": float(result.cost),
        "success": bool(result.success),
        "status": int(result.status),
        "message": str(result.message),
        "nfev": int(result.nfev),
        "optimality": float(result.optimality),
    }


def _json_value(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, quaternion.quaternion):
        return quaternion.as_float_array(value).tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


class BMSProvenanceRecorder:
    """Ordered, checkpointed ledger for replay-relevant frame-fixing state."""

    schema_version = 1

    def __init__(self, output_directory, run_label):
        self.output_directory = Path(output_directory).resolve()
        self.output_directory.mkdir(parents=True, exist_ok=True)
        self.run_label = str(run_label)
        self.events = []
        self.objects = {}
        self.started_utc = datetime.now(timezone.utc).isoformat()

    @property
    def json_path(self):
        return self.output_directory / "bms_provenance_ledger.json"

    @property
    def pickle_path(self):
        return self.output_directory / "bms_provenance_ledger.pkl"

    def record(self, operation_type, objects=None, **details):
        index = len(self.events)
        event = {
            "index": index,
            "operation_type": str(operation_type),
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
            **_json_value(details),
        }
        self.events.append(event)
        if objects:
            for name, value in objects.items():
                self.objects[f"{index}:{name}"] = value
        self.checkpoint()
        return index

    def checkpoint(self):
        ledger = {
            "schema_version": self.schema_version,
            "run_label": self.run_label,
            "started_utc": self.started_utc,
            "events": self.events,
        }
        self._atomic_json(self.json_path, ledger)
        self._atomic_pickle(
            self.pickle_path, {"ledger": ledger, "objects": self.objects}
        )

    @staticmethod
    def _atomic_json(path, value):
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        os.replace(temporary, path)

    @staticmethod
    def _atomic_pickle(path, value):
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("wb") as stream:
            pickle.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, path)
