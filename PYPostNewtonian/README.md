# PYPostNewtonian

Post-Newtonian waveform and evolution code used by `Hybridization.py`.

- `Code/` contains runtime Python modules.  Some are generated from notebooks.
- `PNTerms/` contains symbolic PN source notebooks and polynomial-ratio data.
- `Utilities/` contains helper code used by the notebook generator.

Generated runtime files should be validated by comparing diffs after regeneration.  Notebook outputs
are intentionally not part of the physics source and can be stripped with
`python scripts/strip_notebook_outputs.py` from the repository root.
