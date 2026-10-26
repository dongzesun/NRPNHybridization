# Generated Runtime Code

Runtime modules:

- `PostNewtonian.py` - hand-written API that constructs PN waveforms.
- `PNEvolution.py` - generated PN evolution code.
- `PNWaveformModes.py` - generated strain waveform modes.
- `PNPsiMModes.py` - generated supermomentum/Psi_M modes.
- `EMR22.py` - Schwarzschild flux and binding-energy supplementation through 22PN.
- `EMR22Coefficients.py` and `PNBalance6.py` - generated coefficient tables and legacy balance expressions.

Generator notebooks:

- `PNCodeGen.ipynb`
- `PNCodeGenPsiM.ipynb`

Both notebooks use `GenerateEvolution.ipy` and refresh the supplementary tables
through `scripts/generate_emr22.py` after generating the legacy runtime.

Do not manually reformat generated files unless you are prepared to prove the generated numerical
expressions are unchanged.  See `../../docs/PN_GENERATION.md`.
