#!/usr/bin/env python
"""Strip notebook outputs without changing source cells."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def clean_notebook(path: Path) -> bool:
    nb: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    changed = False
    for cell in nb.get("cells", []):
        if cell.get("cell_type") == "code":
            if cell.get("outputs"):
                cell["outputs"] = []
                changed = True
            if cell.get("execution_count") is not None:
                cell["execution_count"] = None
                changed = True
    if changed:
        path.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", help="Notebook paths; defaults to repository notebooks")
    args = parser.parse_args()

    if args.paths:
        notebooks = [Path(p) if Path(p).is_absolute() else ROOT / p for p in args.paths]
    else:
        notebooks = sorted(ROOT.glob("PYPostNewtonian/**/*.ipynb"))
        notebooks.extend(sorted(ROOT.glob("*.ipynb")))

    changed = [path for path in notebooks if clean_notebook(path)]
    for path in changed:
        print(f"stripped {path.relative_to(ROOT)}")
    print(f"checked {len(notebooks)} notebooks; stripped {len(changed)}")


if __name__ == "__main__":
    main()
