# eccPN

eccPN evolves the relative Cartesian motion of a nonspinning compact binary in
standard harmonic coordinates and computes its time-domain gravitational
radiation. The integration resolves the orbit directly, including its radial
motion. Hereditary forces and radiation use an explicit causal history.

## Installation

From the NRPNHybridization repository root, run

```sh
python -m pip install ./eccPN
```

The package requires Python 3.10 or later, NumPy and SciPy.

## Quick start

This example starts an eccentric binary at periastron in geometric units.
It uses Newtonian conservative dynamics and leading 2.5PN radiation reaction.

```python
import numpy as np
from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve

binary = Binary(0.7, 0.3)
initial = InitialState((25.0, 0.0, 0.0), (0.0, np.sqrt(1.5 / 25.0), 0.0))
trajectory = evolve(binary, initial, None, (0.0, 100.0),
                    accuracy=Accuracy(dynamics="newtonian+2.5PN"))
wave = compute_radiation(trajectory, np.linspace(0.0, 100.0, 101))
h22 = wave.h22
news22 = wave.news22
```

The strain modes represent the asymptotic quantity $R h$. Divide by the
extraction distance to obtain dimensionless strain. The news modes are the
derivatives of the asymptotic strain with respect to retarded time.

## Evolution with a causal past

The default `Accuracy()` uses fixed-history dynamics through 4.5PN with
recoil and three additional 5PN force contributions (failed tail, memory
and reaction squared). The native mass quadrupole and scalar luminosity
include 4PN terms. The mode hierarchy has common 3.5PN amplitude accuracy.
These orders refer to different observables.

```python
import numpy as np
from eccpn import Binary, compute_radiation, evolve
from eccpn.incoming import IncomingOrbit

binary = Binary(0.63, 0.37)
incoming = IncomingOrbit((0.0, 0.18, 0.0), (54.0, 0.0, 0.0), -16000.0,
                         asymptotic_order=3)
initial, past = incoming.prepare(binary)
trajectory = evolve(binary, initial, past, (incoming.cutoff, 1600.0))
wave = compute_radiation(trajectory, np.array([950.0, 1000.0, 1050.0]))
```

The incoming data specify the asymptotic relative velocity, impact parameter
vector and time at which the orbit is matched to its earlier history. Request
radiation on cuts within the evolved trajectory after the matching event.
