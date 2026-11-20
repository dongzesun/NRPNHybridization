# eccPN

`eccPN` computes gravitational waves from two nonspinning compact bodies on
an eccentric orbit. `evolve` follows their motion. `compute_radiation`
returns distance-scaled strain modes and their time derivatives. The
post-Newtonian (PN) model adds corrections in powers of $v/c$ to Newtonian gravity.

Install with `python -m pip install ./eccPN` from the NRPNHybridization root.
Requires Python 3.10 or later, NumPy and SciPy.

## Evolving the orbit

In standard harmonic coordinates, the separation
$\mathbf x=\mathbf y_1-\mathbf y_2$ and relative velocity $\mathbf v$ obey $d\mathbf x/dt=\mathbf v$ and
$d\mathbf v/dt=-Gm\mathbf x/r^3+\mathbf a_{\rm PN}$, where
$m=m_1+m_2$ and $r=|\mathbf x|$. This resolves the rapid motion and burst
of radiation near each closest approach without averaging over an orbit.

The default conservative force goes through 4PN [1]. As waves carry off
energy and angular momentum, the orbit changes. The code includes this
radiation reaction through 4.5PN [2] and the recoil from radiated linear
momentum [3]. Three additional 5PN terms describe the failed tail, orbital
memory and radiation reaction squared. Their starting expressions come
from [4]. We transform them to the same coordinates and time as the orbit.

## Choosing the initial state and past history

Earlier waves scatter off the binary's gravity and act back on the orbit.
This *tail* contributes a history integral to the 4PN acceleration [1]

```math
a^i_{\rm tail,past}(t)=-\frac{8G^2m}{5c^8}x^j(t)
\int_{t_{\rm past}}^t Q^{(7)}_{ij}(s)
\log\!\left[\frac{c(t-s)}{2r(t)}\right]ds .
```

Here $Q_{ij}=\mu(x_i x_j-\delta_{ij}r^2/3)$ is the Newtonian mass quadrupole,
$\mu=m_1m_2/m$, and the superscript counts time derivatives. Two binaries
with the same $\mathbf x$ and $\mathbf v$ today can have different tail forces.

- **Flyby.** Give `IncomingOrbit` the velocity at infinity, the impact-parameter
  vector and a matching time. `prepare` uses an analytic approximation to the
  approach from infinity to compute the position, velocity and accumulated
  recoil at that time, along with the past needed for the radiation.
- **Parabolic approach.** `ParabolicOrbit` uses the pericenter distance and the
  directions of separation and motion at pericenter to construct the incoming orbit.
- **Bound orbit with a finite past.** Give `FiniteOrbitalPast` in `eccpn.history`
  your earlier trajectory, its time boundaries and conserved mass and angular
  momentum coefficients. Supply the initial recoil and center-of-mass integrals.
  The first boundary is $t_{\rm past}$. Earlier emission is absent from the result.

For a physical incoming orbit the lower limit in the tail integral is
$-\infty$. The solver compares one trial step with two half steps. Once they agree
within tolerance, it saves the source as a polynomial on each accepted half
step. Later tail integrals use those polynomials. For a physical incoming
orbit the code integrates over the distant past by parts.

## Extracting the waveforms

At leading order the strain follows the second time derivative of $Q_{ij}$.
The code computes higher mass and current multipoles, adds the propagation
terms, and projects the result onto spin-weighted spherical harmonics [5,6].
For example, the first tail in the radiative quadrupole at retarded time $u$
is [5]

```math
U^{\rm tail}_{ij}(u)=\frac{2GM}{c^3}\int_{t_{\rm past}}^u
Q^{(4)}_{ij}(s)\left[\log\!\left(\frac{u-s}{2b}\right)+\frac{11}{12}\right]ds ,
```

where $M$ is the total gravitating mass and $b$ is the time scale set by `radiative_scale`.
The code returns modes of $H=R(h_+-i h_\times)$.
Divide `wave.h22` by the observer distance $R$ for the dimensionless $(2,2)$
strain. For `wave.news22`, the code differentiates the strain formula
analytically and keeps terms through the requested PN order.

*Nonlinear displacement memory* is the lasting change in strain caused by the energy
the waves carried away. The power radiated in each direction is proportional
to $|N|^2$, where $N=\partial_u H$. Projecting that pattern onto an ordinary
spherical harmonic and integrating gives the change between two times [7]

```math
\Delta H^{\rm mem}_{\ell m}
=\frac{1}{c}\sqrt{\frac{(\ell-2)!}{(\ell+2)!}}
\int_{u_1}^{u_2} (|N|^2)_{\ell m}(u)\,du .
```

An incoming orbit already has accumulated memory at the matching time. The
code adds only terms missing from the published radiative multipoles. This
strain offset is measured at the observer. The orbital-memory force in [4]
acts on the bodies. The modes include terms through 3.5PN, with the mass
quadrupole and total luminosity through 4PN [6].
Memory accuracy here counts powers of $1/c$ at a fixed history. Over a long
inspiral, small errors in the rate can accumulate. A complete 3.5PN
prediction for that accumulated offset would require the secular
radiation reaction through 6PN.

## A flyby in geometric units

Here $G=c=m=1$. The incoming speed is $0.18c$, the impact parameter is
$54Gm/c^2$, and times are in units of $Gm/c^3$.

```python
from eccpn import Binary, compute_radiation, evolve
from eccpn.incoming import IncomingOrbit

binary = Binary(0.63, 0.37)
incoming = IncomingOrbit((0.0, 0.18, 0.0), (54.0, 0.0, 0.0), -16000.0,
                         asymptotic_order=3)
initial, past = incoming.prepare(binary)
trajectory = evolve(binary, initial, past, (incoming.cutoff, 1600.0))
wave = compute_radiation(trajectory, [950.0, 1000.0, 1050.0])
```

The times in the last line sample the burst near closest approach. To
evaluate a different waveform order on the same orbit, pass, for example,
`accuracy=replace(trajectory.accuracy, waveform_pn=3.5)` to
`compute_radiation`, with `replace` from `dataclasses`. A trajectory with
the 3.5PN source history supports 3PN, 3.5PN and the 4PN quadrupole. At
3.5PN the extra memory includes the 4PN quadrupole term. At 4PN that term
is already in the quadrupole.

<details>
<summary>A quick example using just the local force</summary>

This example starts at closest approach. It uses geometric units ($G=c=1$),
Newtonian gravity and the leading 2.5PN radiation-reaction force. This force
depends only on the current position and velocity, so `prehistory=None` is enough.

```python
import numpy as np
from eccpn import Accuracy, Binary, InitialState, compute_radiation, evolve

binary = Binary(0.7, 0.3)
initial = InitialState((25.0, 0.0, 0.0), (0.0, np.sqrt(1.5 / 25.0), 0.0))
trajectory = evolve(binary, initial, None, (0.0, 100.0),
                    accuracy=Accuracy(dynamics="newtonian+2.5PN"))
wave = compute_radiation(trajectory, np.linspace(0.0, 100.0, 101))
h22, news22 = wave.h22, wave.news22
```

</details>

<details>
<summary>Restarting an orbit or changing the output frame</summary>

To restart at the final time `stop`, use
`InitialState(*trajectory.state(stop))` and `trajectory.checkpoint`.
Save with `trajectory.checkpoint.save(path)`. To restore it, import
`AcceptedHistory` from `eccpn.hereditary` and call
`AcceptedHistory.load(path, identity=past.quadrupole.identity)` using the
original source version and evolution settings.

For another Bondi frame, use at least $L+1$ polar and $2L+1$ azimuthal
points. Here $L$ is the largest degree in the native strain and news, the
requested output, or the nonzero supertranslation modes. Check convergence
under angular refinement for boosts and angle-dependent time shifts.

</details>

## References

1. Bernard, Blanchet, Faye & Marchand (2018). Center-of-mass equations of motion at 4PN. Eqs. (2.1)–(2.4) give the local force and Eq. (6.3) the tail. [arXiv:1711.00283](https://arxiv.org/abs/1711.00283v2)
2. Blanchet, Faye, Seraille & Trestini (2026). Radiation reaction in harmonic coordinates. Eqs. (6.6), (6.7) and (6.9)–(6.12) give the 2.5PN, 3.5PN and 4.5PN forces. [arXiv:2601.06743](https://arxiv.org/abs/2601.06743v3)
3. Bini et al. (2026). Radiation reaction on quasi-elliptic and quasi-hyperbolic orbits. Eqs. (2.10a), (5.17) and (5.32b) give the momentum flux and its recoil coupling. [arXiv:2611.12093, v1](https://arxiv.org/abs/2611.12093)
4. Porto, Riva & Yang (2025). Nonlinear gravitational radiation reaction. Sec. 4.3, Eqs. (4.24), (4.25), (4.28) and the supplement give the three 5PN terms in the authors' coordinates. [arXiv:2409.05860](https://arxiv.org/abs/2409.05860v2)
5. Faye, Blanchet & Iyer (2015). Nonlinear multipole interactions. Eqs. (3.4)–(3.10) give the instantaneous terms, tails and memory in the radiative multipoles. [arXiv:1409.3546](https://arxiv.org/abs/1409.3546v2)
6. Blanchet et al. (2023). Gravitational-wave flux and quadrupole modes at 4PN. Sec. II gives the radiative quadrupole and Sec. VI A the general-orbit luminosity. [arXiv:2304.11186](https://arxiv.org/abs/2304.11186v4)
7. Favata (2009). Post-Newtonian corrections to gravitational-wave memory. Eqs. (2.26) and (2.31) relate angular energy flux to displacement memory. [arXiv:0812.0069](https://arxiv.org/abs/0812.0069v2)

<details open>
<summary>Where each term comes from, including coefficients derived for eccPN</summary>

### Orbit and initial state

| Quantity | Formula used in the code |
| --- | --- |
| Relative acceleration | The local conservative coefficients and their logarithmic length scales are [1], Eqs. (2.1)–(2.4). The 4PN tail, including the two local terms that accompany the integral above, is the retarded relative acceleration obtained from [1], Eq. (6.3). The 2.5PN, 3.5PN and 4.5PN reaction terms use [2], Eqs. (6.6), (6.7), (6.9)–(6.12). [Local force](src/eccpn/local_force.py) · [Tail force](src/eccpn/tails.py) |
| Recoil | The radiated momentum and center-of-mass integral obey $\dot\Pi=\mathcal F_P$ and $\dot\Gamma=\Pi+\mathcal F_G$ from [2], Eqs. (6.1)–(6.4), with fluxes from [3], Eqs. (2.10a,b), (5.17). In [2], Eq. (6.10), the 4.5PN force uses $\Pi_{\rm lead}+\Pi_{\rm next}+\mathcal F_G$. The extra coupling in [3], Eq. (5.32b), uses only $\Pi_{\rm lead}$ and can be enhanced on bound orbits. [Recoil](src/eccpn/recoil.py) |
| Body positions and conserved quantities | The positions use [1], Eqs. (B4)–(B6), and [2], Eqs. (6.2)–(6.4). The waves use the conserved mass and angular momentum of [1], Sec. III, and [5], Sec. II. For a flyby we calculate these at infinity from the input speed and impact parameter. A finite history supplies its own coefficients. The native output reports the enhanced recoil correction to the quadrupole separately. Pointwise balance checks use the Schott terms in [1], Eqs. (6.6)–(6.11), and [2], Eqs. (6.13)–(6.15). [Worldlines](src/eccpn/worldlines.py) · [Recoil in the quadrupole](src/eccpn/recoil_radiation.py) |
| Three 5PN forces | We start from `aFT`, `aM` and `aRR2` in the supplement to [4]. We derive the change of position and clock to standard harmonic coordinates, including its action on the lower-order force. Evolving the 2.5PN force also changes the trajectory at second order. The analytic checks use the static-spin metric, the quadrupole-squared metric and the independent scattering results below. [Harmonic coefficients](src/eccpn/generated/radiative_reaction_5pn.py) · [Henry & Larrouturou, Eqs. (3.9)–(3.12)](https://arxiv.org/abs/2307.05860v3) · [Blanchet, Eq. (4.13) and Table 2](https://arxiv.org/abs/gr-qc/9710037v1) · [Porto & Riva, supplement](https://arxiv.org/abs/2604.09545v3) · [Brunello, Sec. V](https://arxiv.org/abs/2610.42057v1) |
| Incoming history and matching | We expand the harmonic force [1,2] on the incoming branch. The hyperbolic expansion uses inverse time and logarithms. The parabolic expansion starts from the Newtonian parabolic orbit. Earlier momentum flux and memory supply the integration constants at the join. To match a specified state, the shooting routine varies incoming data and evolves forward. The expansion, matching and quadrature are calculations in this code. [Incoming trajectories](src/eccpn/incoming.py) · [Incoming radiation](src/eccpn/incoming_waveform.py) |

### The mass and current moments

Here $I_\ell$ and $J_\ell$ mean the source mass and current tensors of rank
$\ell$. The common 3.5PN waveform reaches $I_9$ and $J_8$. A term with
an explicit $c^{-k}$ enters the strain at $(\ell-2+k)/2$ PN for mass or
$(\ell-1+k)/2$ PN for current.

| Quantity | Published input and derivation |
| --- | --- |
| Leading moments and published eccentric coefficients | The Newtonian moments at every rank follow [Blanchet et al., Eqs. (5.21)–(5.22)](https://arxiv.org/abs/0802.1249v3). The general-orbit $I_2$ through 3.5PN is in [Faye et al., Eqs. (A1)–(A3)](https://arxiv.org/abs/1204.1043v2). $I_3$ uses [5], Eqs. (4.9)–(4.10), and $J_2$ uses [Henry et al., Eqs. (3.6)–(3.7)](https://arxiv.org/abs/2105.10876v2). The printed higher-rank coefficients come from [Mishra et al., Eqs. (3.23), (3.26)](https://arxiv.org/abs/1501.07096v2). [Leading moments](src/eccpn/leading_modes.py) · [Higher moments](src/eccpn/higher_moments.py) |
| Higher-rank coefficients derived here | At common waveform order 3.5PN, the odd terms in $I_4$ and $J_3$ use the retarded source and the harmonic center-of-mass shift. $I_5$ uses the general 2PN mass integral. $J_4$ uses the compact and quadratic field terms. The circular $(4,1)$ and $(4,3)$ modes fix its one remaining cubic coefficient. The 1PN part of $J_6$ follows from the current integral. The same mass integral at 1PN also gives the correction to $I_7$. [Coefficients](src/eccpn/generated/higher_multipoles.py) · [Mass and current integrals, Eqs. (5.15)–(5.16)](https://arxiv.org/abs/gr-qc/9801101v1) · [Odd mass term, Eq. (4.12)](https://arxiv.org/abs/gr-qc/9603048v2) · [2PN mass integral, Eqs. (3.5), (3.28), (3.69b)](https://arxiv.org/abs/gr-qc/9501029) · [Field source, Eqs. (3), (5b), (54)–(55)](https://arxiv.org/abs/gr-qc/0105098v2) · [Henry, circular modes](https://arxiv.org/abs/2210.15602v2) |
| General-orbit 4PN source quadrupole | The source integral, nonlocal part and finite renormalization come from the three papers below. We include the infrared commutator of [6], Eq. (3.24), once. We reconstruct the remaining local tensor coefficients at equal harmonic ultraviolet scales from the general-orbit scalar `FluxCanLoc4PN` in the supplement to [6], together with the independently calculated terms linear in $G$ and the known logarithms. The circular tensor supplies a check. [Local coefficients](src/eccpn/generated/local_quadrupole_4pn.py) · [Source integral](https://arxiv.org/abs/2003.13672v2) · [Nonlocal part, Eq. (2.14)](https://arxiv.org/abs/2110.02240v3) · [Finite renormalization, Eq. (6.10)](https://arxiv.org/abs/2110.02243v3) |

### Propagation and observer quantities

| Quantity | Formula used in the code |
| --- | --- |
| Source to canonical moments | The gauge moments $W,X,Y,Z$ are from [2], Appendix D, Eqs. (D4a–d,g,i,k). The lower-order map is [5], Eqs. (3.11a–e). For the 4PN quadrupole we use the dimensionally regularized relation [6], Eq. (2.39). [Quadrupole map](src/eccpn/canonical.py) · [Higher-rank maps](src/eccpn/mpm.py) |
| Canonical to radiative moments, strain and news | Instantaneous products, ordinary tails, tails of tails and the explicit memory use [5], Eqs. (3.4)–(3.10). The 4PN quadrupole also includes the cubic terms in [6], Eq. (2.7), with their ordered integrals. Each time derivative uses the equations of motion truncated at the order needed for that moment. We then project onto spin-weight $-2$ harmonics. The time argument is the matched retarded time in [5], Eq. (3.1). Its scale $b$ is a time, distinct from the harmonic coordinate lengths and the infrared length $r_0$. [Moment derivatives](src/eccpn/waveform_jets.py) · [Radiative quadrupole](src/eccpn/radiative_quadrupole.py) · [Tail kernels](src/eccpn/radiation_kernels.py) · [Cubic tail of memory](https://arxiv.org/abs/2301.09395v2) · [Kidder, projection Eqs. (13)–(19)](https://arxiv.org/abs/0710.0614v1) |
| Energy and angular momentum at infinity | The scalar luminosity has the three parts in [6], Eqs. (6.1)–(6.5). The code evaluates the nonlocal Eq. (6.4) as printed and also computes a version using the actual retarded orbital tail and its derivative. The default sums this retarded version, the published local scalar and nonlinear propagation. Both nonlocal values are available in `wave.energy_flux_parts` for native-frame 4PN output. The code also contracts the mode news for power and the radiative moments and their derivatives for angular momentum. The scalar supplied one constraint on the reconstructed source quadrupole above, so agreement between scalar and mode power tests their assembly. [Luminosity and torque](src/eccpn/luminosity.py) · [Thorne, Sec. IV, especially Eq. (4.23′)](https://doi.org/10.1103/RevModPhys.52.299) |
| Memory | The angular products in the flux integral above use [7], Eqs. (2.26), (2.31). The boundary terms relating this integral to the radiative mass moments follow [5], Appendix C. We use the same past for both and subtract the overlapping orders. Products of news through the 3.5PN summed order generate memory modes up to degree 11. Ordinary scattering memory is already in the difference of the two endpoint strains. The current moments carry the magnetic strain associated with spin memory through [5], Eqs. (3.5b,d). [Memory accumulation](src/eccpn/memory.py) · [Nichols, Sec. IV D](https://arxiv.org/abs/1702.03300) |
| Observer frame | We use $u_g=cu$ and shear $\sigma=\overline H/2$. A boost, rotation or supertranslation changes both the observer direction and the time sampled along it. The incoming mass aspect uses the velocities at infinity. [Frame transformation](src/eccpn/frame.py) · [Mitman et al., Eqs. (7), (9), (13), (29)](https://arxiv.org/abs/2208.04356) · [Riva et al., Eqs. (3.9)–(3.13)](https://arxiv.org/abs/2302.09065) |

The code implements the finite-start tail by integrating the stored
$Q^{(7)}$ directly on each polynomial panel. Differentiating that integral
also differentiates its moving endpoint and its scale. Saving accepted
steps, checking restart identity, evaluating powers of $1/c$ in floating-point
arithmetic in the chosen units, choosing compatible stored moments for a requested PN order and
checking the angular grid are numerical parts of the implementation.
[Panel integrals](src/eccpn/history.py) · [Tail derivative](src/eccpn/tail_derivative.py) · [Accepted history](src/eccpn/hereditary.py) · [Radiation orders](src/eccpn/radiation.py)

For comparisons in modified harmonic coordinates we use the 3PN shift
of [Arun et al., Eqs. (4.3)–(4.10)](https://arxiv.org/abs/0711.0302v2).
We check the tails and memory against independent eccentric calculations
and the published circular limit.
[Eccentric tails](https://arxiv.org/abs/1904.11814v3) · [Eccentric memory](https://arxiv.org/abs/1906.06263v2) · [Circular memory](https://arxiv.org/abs/2410.23950v2)

</details>
