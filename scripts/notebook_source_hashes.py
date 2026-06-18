#!/usr/bin/env python
"""Print hashes of notebook source cells, ignoring outputs and execution counts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def source_hash(path: Path) -> str:
    nb = json.loads(path.read_text(encoding="utf-8"))
    source = "".join("".join(cell.get("source", [])) for cell in nb.get("cells", []))
    return hashlib.sha256(source.encode("utf-8")).hexdigest()


def main() -> None:
    notebooks = sorted(ROOT.glob("PYPostNewtonian/**/*.ipynb"))
    notebooks.extend(sorted(ROOT.glob("*.ipynb")))
    for path in notebooks:
        print(f"{source_hash(path)}  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
