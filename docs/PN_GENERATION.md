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

## Schwarzschild supplementation through 22PN

`Code/EMR22.py` evaluates the additional nonspinning test-particle terms without
expanding the spin/precession/mode generators to unsupported orders. Both
`PNCodeGen` notebooks share `Code/GenerateEvolution.ipy` for the evolution emitter,
so the Psi_M generator retains the flux API. They refresh the derived balance
helper and coefficient tables after writing the legacy runtime.

Regenerate just the new numerical data with:

```bash
python scripts/generate_emr22.py
```

This requires mpmath, and emits `Code/EMR22Coefficients.py`, `Code/PNBalance6.py`,
and the decimal reference table. The 6PN balance helper is extracted from the
existing generated TaylorT1 expressions; the original source notebooks remain
the authority for these expressions. The exact Fujita source and its provenance
are retained under `PNTerms/EMR22`. The supplement's dimensionless series is
multiplied by `(32/5)*nu**2*v**10`; the README's single power of mu/M is not used
as the physical total-flux normalization.

The original runtime expressions through 6PN are unchanged by this extension.
Only powers beyond relative 6PN are added. The binding energy receives the same
orders from the Schwarzschild circular-orbit expression. Horizon absorption,
spin/precession and waveform-amplitude orders remain unchanged. These are
EMR-supplemented comparable-mass models, not complete 22PN dynamics.

For orders above 6, the TaylorT1 balance equations are integrated using velocity
as the independent variable. This is a change of integration variable, not a
re-expansion of the TaylorT1 ratio. Integration stops at the energy turning point;
the returned waveform uses the historical Newtonian time-grid estimate, bounded
below by the shortest orbital timestep on the integrated branch. This avoids
vanishing output steps at a Newtonian merger time that differs from the PN endpoint.

Example:

```python
PNWaveform(q, omega, chi1, chi2,
           PNEvolutionOrder=22, PNWaveformModeOrder=6)
```

The CLI equivalents are `--PNEvolutionOrder 22 --PNWaveformModeOrder 6`.
`--PNOrder` retains its legacy meaning and defaults. Direct flux-on-state-array
calls use `PNEvolutionOrder` as the flux truncation order without evolving.

Numerical checks: `python tests/test_emr22.py` with this repository on PYTHONPATH.
