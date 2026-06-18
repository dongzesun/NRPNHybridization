# NRPNHybridization

Tools for aligning and hybridizing numerical-relativity waveforms with post-Newtonian waveforms.

The main driver is `Hybridization.py`.  The post-Newtonian waveform implementation lives in
`PYPostNewtonian/Code`; those Python files are generated from notebooks and should be treated as
physics outputs, not hand-edited implementation code.

## Repository Layout

- `Hybridization.py` - command-line driver for NR/CCE to PN alignment and hybridization.
- `PNBMS.py` - BMS-frame mapping helpers used by the driver.
- `PYPostNewtonian/Code/` - generated PN evolution and waveform modules used at runtime.
- `PYPostNewtonian/PNTerms/` - symbolic PN terms and source notebooks used by the code generator.
- `PYPostNewtonian/Utilities/` - notebook/code-generation helpers.
- `docs/` - notes on structure and the generation workflow.
- `scripts/` - maintenance helpers that avoid editing physics content by hand.

## Running Hybridization

```bash
python Hybridization.py \
  --t -24000.0 \
  --SimDir /path/to/simulation/LevN \
  --CCEDir /path/to/cce/LevN \
  --OutName Out12PNS24000L5.npz \
  --nOrbits 5.0 \
  --ABDFormat CCE
```

Useful arguments:

- `--t`: end time of the matching window.
- `--SimDir`: directory containing horizon/simulation data.
- `--CCEDir`: directory containing CCE waveform data.
- `--OutName`: checkpoint/result `.npz` output.  If this file exists, the script resumes from it.
- `--length`: matching-window length, ignored when `--nOrbits` is given.
- `--nOrbits`: matching-window length in NR orbits.
- `--truncate t1 t2`: keep only CCE data between `t1` and `t2` to reduce memory use.
- `--ABDFormat`: format label passed to `scri.SpEC.create_abd_from_h5`; the historical default is `RPDMB`, while `CCE` is required by some newer `scri` environments.

## Physics-Code Policy

Do not hand-edit generated files in `PYPostNewtonian/Code` unless you are intentionally changing the
physics and can regenerate/validate the outputs.  Prefer changing the notebooks in
`PYPostNewtonian/PNTerms` or `PYPostNewtonian/Code/PNCodeGen*.ipynb`, then regenerate the output files.

See:

- `docs/STRUCTURE.md`
- `docs/PN_GENERATION.md`
