# PN Code Generation

The post-Newtonian implementation uses notebooks as symbolic sources and generated Python files as
runtime code.

## Generated Runtime Files

- `PYPostNewtonian/Code/PNEvolution.py`
- `PYPostNewtonian/Code/PNWaveformModes.py`
- `PYPostNewtonian/Code/PNPsiMModes.py`

`PYPostNewtonian/Code/PostNewtonian.py` is hand-written orchestration code that imports the generated
modules.

## Generator Notebooks

- `PYPostNewtonian/Code/PNCodeGen.ipynb` generates `PNEvolution.py` and `PNWaveformModes.py`.
- `PYPostNewtonian/Code/PNCodeGenPsiM.ipynb` generates `PNEvolution.py` and `PNPsiMModes.py`.

Both notebooks execute source notebooks under `PYPostNewtonian/PNTerms`.  The symbolic inputs are
messy but should be treated as the source of truth for physics expressions.

## Recommended Workflow

1. Record source hashes:

   ```bash
   python scripts/notebook_source_hashes.py
   ```

2. Edit only the relevant source notebook(s).

3. Regenerate generated runtime code:

   ```bash
   python scripts/regenerate_pn_code.py
   ```

4. Inspect generated diffs carefully:

   ```bash
   git diff -- PYPostNewtonian/Code/PNEvolution.py               PYPostNewtonian/Code/PNWaveformModes.py               PYPostNewtonian/Code/PNPsiMModes.py
   ```

5. Strip notebook outputs before committing:

   ```bash
   python scripts/strip_notebook_outputs.py
   ```

## Compatibility Note

The notebooks use historical IPython idioms such as `execfile` and notebook execution helpers in
`PYPostNewtonian/Utilities/ExecNotebook.ipy`.  If regeneration fails in a modern Python environment,
fix the execution harness first and verify that generated `.py` output is unchanged before changing
any physics expressions.
