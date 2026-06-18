# Repository Structure

This repository intentionally keeps the historical top-level module layout intact so existing scripts
that import `PNBMS`, `Hybridization`, or `PYPostNewtonian` continue to work.

## Runtime Path

`Hybridization.py` imports:

1. `scri`, `numpy`, `scipy`, `quaternion`, and related scientific packages.
2. `PNBMS.py` for BMS-frame mapping.
3. `PYPostNewtonian.Code.PostNewtonian`, which calls the generated PN modules.

The runtime PN modules are:

- `PYPostNewtonian/Code/PostNewtonian.py`
- `PYPostNewtonian/Code/PNEvolution.py`
- `PYPostNewtonian/Code/PNWaveformModes.py`
- `PYPostNewtonian/Code/PNPsiMModes.py`

## Source-Generation Path

The PN source notebooks are in `PYPostNewtonian/PNTerms`.  The generator notebooks in
`PYPostNewtonian/Code` assemble those symbolic terms into the generated Python files.

Important generator inputs:

- `PYPostNewtonian/PNTerms/Variables_Q.ipynb`
- `PYPostNewtonian/PNTerms/OrbitalEvolution.ipynb`
- `PYPostNewtonian/PNTerms/AngularMomentum.ipynb`
- `PYPostNewtonian/PNTerms/WaveformModes.ipynb`
- `PYPostNewtonian/PNTerms/PsiMModes.ipynb`
- `PYPostNewtonian/PNTerms/PolynomialRatios/*.dat`

## Maintenance Rules

- Keep generated `.py` files committed so production runs do not depend on notebook execution.
- Strip notebook outputs before committing; outputs are environment-specific and not physics content.
- Use `scripts/notebook_source_hashes.py` before and after cleanup to verify source cells did not change.
- Use `scripts/regenerate_pn_code.py --dry-run` to see which notebooks and outputs are involved.

## Runtime Compatibility

`Hybridization.py` exposes `--ABDFormat` so the `scri.SpEC.create_abd_from_h5` format label can be selected without editing source code.  The default remains the historical `RPDMB`; use `--ABDFormat CCE` on environments where current `scri` expects CCE block-mode files.
