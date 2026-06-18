# PN Term Source Notebooks

These notebooks define symbolic post-Newtonian terms used by the generator notebooks in
`../Code`.

Typical dependency flow:

1. `Variables_Q.ipynb` defines the basic symbolic variables.
2. `BindingEnergy.ipynb`, `Flux.ipynb`, `Precession.ipynb`, and `EnergyAbsorption.ipynb` define source terms.
3. `OrbitalEvolution.ipynb` combines source terms and polynomial-ratio data.
4. `WaveformModes.ipynb` and `PsiMModes.ipynb` define mode expressions.
5. `../Code/PNCodeGen*.ipynb` writes generated Python modules.

Keep outputs stripped; code and markdown cells are the source content.
