#!/usr/bin/env python
"""Run the PN code-generation notebooks from a controlled entry point."""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "PYPostNewtonian" / "Code"
OUTPUTS = [
    CODE_DIR / "PNEvolution.py",
    CODE_DIR / "PNWaveformModes.py",
    CODE_DIR / "PNPsiMModes.py",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "missing"


def run_notebook(notebook: str) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = f"{ROOT}:{ROOT / 'PYPostNewtonian'}:{env.get('PYTHONPATH', '')}"
    cmd = [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        "--inplace",
        notebook,
    ]
    subprocess.run(cmd, cwd=CODE_DIR, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="Print actions without executing notebooks")
    parser.add_argument(
        "--psi-m-only",
        action="store_true",
        help="Only run PNCodeGenPsiM.ipynb instead of both generator notebooks",
    )
    args = parser.parse_args()

    notebooks = ["PNCodeGenPsiM.ipynb"] if args.psi_m_only else ["PNCodeGen.ipynb", "PNCodeGenPsiM.ipynb"]

    print("generator notebooks:")
    for notebook in notebooks:
        print(f"  {CODE_DIR / notebook}")

    print("generated outputs before:")
    before = {path: sha256(path) for path in OUTPUTS}
    for path, digest in before.items():
        print(f"  {digest}  {path.relative_to(ROOT)}")

    if args.dry_run:
        return

    for notebook in notebooks:
        run_notebook(notebook)

    print("generated outputs after:")
    for path in OUTPUTS:
        digest = sha256(path)
        changed = "changed" if digest != before[path] else "unchanged"
        print(f"  {digest}  {path.relative_to(ROOT)}  {changed}")


if __name__ == "__main__":
    main()
