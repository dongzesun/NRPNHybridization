"""Fixed-history radiative assembly and the higher formal memory source."""

from dataclasses import dataclass, replace
from functools import cached_property, partial

import numpy as np

from .compact_moments import CompactMomentHistory, pack_stf
from .generated import mpm_tables
from .hereditary import _interpolate
from .history import PanelHistory, SourcePanel
from .incoming_waveform import IncomingMomentHistory, IncomingWaveform
from .leading_modes import leading_mode_coefficients
from .mpm import _contraction_equation, _LEVI_CIVITA, instantaneous_corrections
from .multipoles import known_local_quadrupole_grade_eight, nonlocal_quadrupole_grade_eight
from .projection import project_stf, stf_tensor
from .radiative_quadrupole import mpm_quadrupole_interactions
from .radiation_kernels import ordinary_tail, tail_of_tail
from .waveform_jets import canonical_waveform_jets


@dataclass(frozen=True)
class CanonicalHistory:
    """Canonical jets over one finite history and its conserved mass grades."""

    moments: dict
    mass_grades: dict
    dipole_grades: dict
    harmonic_scale: float
    infrared_scale: float
    maximum_absolute_grade: int = 6
    incoming: IncomingWaveform | None = None

    def covered_moments(self):
        if self.incoming is None:
            return self.moments
        return {name: {grade: IncomingMomentHistory(past, self.incoming, name, grade)
                       for grade, past in grades.items()} for name, grades in self.moments.items()
                if name != "orbital"}

    @cached_property
    def quadrupole_for_cubic_terms(self):
        finite = self.moments["M2"][0]
        if self.incoming is None:
            return finite
        from .asymptotic_waveform import AsymptoticQuadrupole, MatchedQuadrupoleHistory, SampledWaveformQuadrupole
        from .incoming import IncomingOrbit

        orbit = self.incoming.orbit
        is_hyperbolic = isinstance(orbit,IncomingOrbit)
        asymptote = (AsymptoticQuadrupole.hyperbolic_newtonian if is_hyperbolic
                     else AsymptoticQuadrupole.parabolic_newtonian)
        past = asymptote(self.incoming.binary,orbit)
        complete_order = 3 if is_hyperbolic else 2
        if orbit.asymptotic_order==complete_order:
            past = SampledWaveformQuadrupole(self.incoming.binary,orbit,past,self.incoming.order)
        return MatchedQuadrupoleHistory(finite, past)

    def extend(self, other):
        if (self.harmonic_scale != other.harmonic_scale or self.infrared_scale != other.infrared_scale
                or self.maximum_absolute_grade != other.maximum_absolute_grade
                or self.mass_grades != other.mass_grades
                or self.dipole_grades.keys() != other.dipole_grades.keys()
                or any(not np.array_equal(charge, other.dipole_grades[grade])
                       for grade, charge in self.dipole_grades.items())
                or self.moments.keys() != other.moments.keys()
                or any(grades.keys() != other.moments[name].keys() for name, grades in self.moments.items())):
            raise ValueError("Canonical history extension requires the same moments, charges and scales.")
        moments = {
            name: {
                grade: (past.extend(other.moments[name][grade]) if isinstance(past,CompactMomentHistory)
                        else past.append(*other.moments[name][grade].panels))
                for grade,past in grades.items()
            }
            for name, grades in self.moments.items()
        }
        return CanonicalHistory(moments, self.mass_grades, self.dipole_grades, self.harmonic_scale,
                                self.infrared_scale, self.maximum_absolute_grade, self.incoming)

    def to_record(self):
        panels = self.moments["M2"][0].panels
        return {
            "cuts": [panels[0].start, *(panel.end for panel in panels)],
            "masses": self.mass_grades,
            "dipoles": {grade: np.asarray(charge).tolist() for grade, charge in self.dipole_grades.items()},
            "harmonic": self.harmonic_scale,
            "infrared": self.infrared_scale,
            "grade": self.maximum_absolute_grade,
            "incoming": None if self.incoming is None else self.incoming.to_record(),
            "moments": {
                name: {
                    grade: [panel.coefficients.tolist() for panel in (
                        past.packed.panels if isinstance(past,CompactMomentHistory) else past.panels
                    )] for grade,past in grades.items()
                }
                for name, grades in self.moments.items()
            },
            "packed_ranks": {name: next(iter(grades.values())).rank for name,grades in self.moments.items()
                             if isinstance(next(iter(grades.values())),CompactMomentHistory)},
        }

    @classmethod
    def from_record(cls, record):
        cuts = record["cuts"]

        def restore(name, grade, coefficients):
            if len(coefficients) != len(cuts) - 1 or not coefficients:
                raise ValueError("Every saved canonical jet must cover all shared cuts.")
            panels = tuple(SourcePanel(a, b, data) for a, b, data in zip(cuts, cuts[1:], coefficients))
            finite = PanelHistory(cuts[0],panels[0].coefficients.shape[1:],f"{name}-grade-{grade}",panels)
            packed_rank = record.get("packed_ranks",{}).get(name)
            return CompactMomentHistory(finite,packed_rank) if packed_rank is not None else finite

        incoming = record.get("incoming")
        return cls(
            {name: {int(g): restore(name, g, rows) for g, rows in grades.items()}
             for name, grades in record["moments"].items()},
            {int(grade): value for grade, value in record["masses"].items()},
            {int(grade): np.asarray(value) for grade, value in record["dipoles"].items()},
            record["harmonic"], record["infrared"], record["grade"],
            None if incoming is None else IncomingWaveform.from_record(incoming),
        )


def build_canonical_history(
    binary, times, state, mass_grades, dipole_grades, harmonic_scale, infrared_scale,
    maximum_absolute_grade=6, incoming=None,
):
    """Sample an explicit covered SH state into quadratic canonical panels."""
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all() or not (np.diff(times) > 0).all():
        raise ValueError("Supply two or more strictly increasing finite history boundaries.")
    expected_mass = {0, 2} if maximum_absolute_grade == 6 else {0, 2, 4}
    expected_dipole = {0} if maximum_absolute_grade == 6 else {0, 2}
    if (maximum_absolute_grade not in (6, 7) or not expected_mass <= mass_grades.keys()
            or not expected_dipole <= dipole_grades.keys() or mass_grades[0] != binary.mass or any(
            grade not in (0, 2, 4) or not np.isfinite(value) for grade, value in mass_grades.items()
    )):
        raise ValueError("Supply all finite conserved mass and dipole grades required by the waveform budget.")

    def sample(time):
        position, velocity = state(time)
        jets = canonical_waveform_jets(
            binary, position, velocity, harmonic_scale, infrared_scale, dipole_grades,
            maximum_absolute_grade,
        )
        for name,grades in jets.items():
            if _stored_rank(name) is not None:
                jets[name] = {grade:pack_stf(jet,_stored_rank(name)) for grade,jet in grades.items()}
        jets["orbital"] = {0: np.stack((position, velocity))}
        return jets

    output = {}
    initial = sample(times[0])
    for start, stop in zip(times[:-1], times[1:]):
        middle, final = sample((start + stop) / 2), sample(stop)
        for name, inventory in initial.items():
            panels = output.setdefault(name, {})
            for grade, first in inventory.items():
                coefficients = _interpolate((first, middle[name][grade], final[name][grade]))
                panels.setdefault(grade, []).append(SourcePanel(start, stop, coefficients))
        initial = final
    def history(name,grade,panels):
        finite = PanelHistory(
            times[0],panels[0].coefficients.shape[1:],f"{name}-grade-{grade}",tuple(panels),
        )
        rank = _stored_rank(name)
        return finite if rank is None else CompactMomentHistory(finite,rank)

    return CanonicalHistory({
        name:{grade:history(name,grade,panels) for grade,panels in data.items()}
        for name,data in output.items()
    },dict(mass_grades),dict(dipole_grades),harmonic_scale,infrared_scale,maximum_absolute_grade,incoming)


def _stored_rank(name):
    if name.startswith(("M","S")) and int(name[1:])>=3:
        return int(name[1:])
    return None


def _add(pieces, grade, value):
    pieces[grade] = pieces.get(grade, 0) + value


@dataclass(frozen=True)
class _FormalMemoryProduct:
    cuts: np.ndarray
    coefficients: np.ndarray
    prefix: np.ndarray

    @classmethod
    def from_moments(cls,first,second,tensors,equation):
        orders = [(symbol,order) for symbol,order,_ in tensors if symbol!="epsilon"]
        first_panels,second_panels = first.panels,second.panels
        cuts = np.array([first.start,*(panel.end for panel in first_panels)])
        other_cuts = np.array([second.start,*(panel.end for panel in second_panels)])
        if not np.array_equal(cuts,other_cuts):
            raise ValueError("Formal MPM memory sources must share the same accepted panels.")
        factors,output = equation.split("->")
        labels = ["z"+label if symbol!="epsilon" else label
                  for label,(symbol,_,_) in zip(factors.split(","),tensors)]
        contraction = ",".join(labels)+"->z"+output

        def product(left,right):
            variables = iter((left,right))
            arguments = [_LEVI_CIVITA if symbol=="epsilon" else next(variables) for symbol,_,_ in tensors]
            return np.einsum(contraction,*arguments,optimize=True)

        left = np.stack([panel.coefficients[:,orders[0][1]] for panel in first_panels])
        right = np.stack([panel.coefficients[:,orders[1][1]] for panel in second_panels])
        shape = product(left[:,0],right[:,0]).shape
        polynomial = np.zeros((shape[0],left.shape[1]+right.shape[1]-1,*shape[1:]))
        for i in range(left.shape[1]):
            for j in range(right.shape[1]):
                polynomial[:,i+j] += product(left[:,i],right[:,j])
        degree = np.arange(1,polynomial.shape[1]+1)
        integral = np.einsum("d,zd...->z...",1/degree,polynomial)
        integral *= np.diff(cuts).reshape((-1,)+(1,)*(integral.ndim-1))
        incoming = np.zeros(shape[1:])
        if isinstance(first,IncomingMomentHistory):
            if not isinstance(second,IncomingMomentHistory) or first.past is not second.past:
                raise ValueError("Formal memory products must share the incoming quadrature grid.")
            a = first.past.sampled[first.name][first.grade][:,orders[0][1]]
            b = first.past.sampled[second.name][second.grade][:,orders[1][1]]
            incoming = np.tensordot(first.past.grid[1],product(a,b),axes=1)
        prefix = np.concatenate((incoming[None,...],incoming+np.cumsum(integral,axis=0)))
        return cls(cuts,polynomial,prefix)

    def value(self,time):
        index = np.searchsorted(self.cuts,time,side="left")-1
        if index<0:
            return self.prefix[0]
        width = self.cuts[index+1]-self.cuts[index]
        fraction = (time-self.cuts[index])/width
        order = np.arange(1,self.coefficients.shape[1]+1)
        return self.prefix[index]+width*np.tensordot(
            fraction**order/order,self.coefficients[index],axes=1)


def _panel_products(first,second,tensors,equation,time,cache):
    key = (id(getattr(first,"finite",first)),id(getattr(second,"finite",second)),tensors,equation)
    if key not in cache:
        cache[key] = _FormalMemoryProduct.from_moments(first,second,tensors,equation)
    return cache[key].value(time)


def _memory_value(binary,target,canonical,time,maximum_relative_grade,cache):
    result = {}
    rank = int(target[1])
    for term_grade, numerator, denominator, tensors in mpm_tables.MEMORY.get(target, ()):
        if term_grade > maximum_relative_grade:
            continue
        names = [symbol for symbol, _, _ in tensors if symbol != "epsilon"]
        equation = _contraction_equation(tensors, rank)
        scale = binary.G * numerator / (binary.c**term_grade * denominator)
        for grade_a, first in canonical[names[0]].items():
            for grade_b, second in canonical[names[1]].items():
                grade = term_grade + grade_a + grade_b
                if grade > maximum_relative_grade:
                    continue
                value = _panel_products(first,second,tensors,equation,time,cache)
                _add(result, grade, scale * value)
    return {grade: stf_tensor(value, rank) for grade, value in result.items()}


def finite_radiation_available(
    binary, history, time, position, velocity, radiative_scale, *,
    news_only=False, sampled_jets=None, news_endpoint=None, log_history=None, prefix_cache=None,
):
    """Graded native modes, news and formal memory on the specified finite past."""
    from .mpm import formal_memory_rates

    if prefix_cache is None:
        prefix_cache = {}
    memory_prefix = prefix_cache.setdefault(("formal_memory",id(history)),{})
    maximum = history.maximum_absolute_grade
    canonical = history.covered_moments()
    current = sampled_jets
    if current is None:
        current = canonical_waveform_jets(
            binary, position, velocity, history.harmonic_scale, history.infrared_scale,
            history.dipole_grades, maximum,
        )
    derivatives = (1,) if news_only else (0, 1)
    result = ({}, {})
    tensors = {}
    for kind, targets in (("mass", "u"), ("current", "v")):
        explicit_rank = 7 if maximum == 7 and kind == "mass" else 6
        for rank in range(2, explicit_rank + 1):
            offset = rank - 2 + (kind == "current")
            if offset > maximum:
                continue
            name = ("M" if kind == "mass" else "S") + str(rank)
            target = targets + str(rank)
            budget = maximum - offset
            pieces = ({}, {})
            for grade, jet in current[name].items():
                if grade <= budget:
                    pieces[0][grade], pieces[1][grade] = jet[rank:rank + 2]
            for derivative in derivatives:
                additions = instantaneous_corrections(binary, target, current, derivative, maximum)
                for grade, tensor in additions.items():
                    _add(pieces[derivative], grade, tensor)
            if kind == "mass":
                rates = formal_memory_rates(binary, target, current, maximum)
                for grade, rate in rates.items():
                    if grade <= budget:
                        _add(pieces[1], grade, rate)
                if not news_only:
                    for grade, memory in _memory_value(binary,target,canonical,time,budget,memory_prefix).items():
                        _add(pieces[0], grade, memory)
            for grade, past in canonical[name].items():
                for mass_grade, mass in history.mass_grades.items():
                    additional = grade + mass_grade + 3
                    if additional <= budget:
                        endpoint = {} if isinstance(news_endpoint, dict) else news_endpoint
                        tail = ordinary_tail(binary, past, time, rank, kind, mass, radiative_scale,
                                             endpoint=endpoint,
                                             split_jet=None if log_history is None else partial(log_history.jet,name,grade,parts=2))
                        if isinstance(endpoint, dict):
                            _project_endpoint(binary, kind, rank, additional + offset, endpoint, news_endpoint)
                        for derivative in derivatives:
                            _add(pieces[derivative], additional, tail[derivative])
            if budget >= 6:
                endpoint = {} if isinstance(news_endpoint, dict) else news_endpoint
                tail = tail_of_tail(
                    binary, canonical[name][0], time, rank, kind, history.mass_grades[0],
                    radiative_scale, history.infrared_scale, endpoint=endpoint,
                    split_jet=None if log_history is None else partial(log_history.jet,name,0,parts=2),
                )
                if isinstance(endpoint, dict):
                    _project_endpoint(binary, kind, rank, 6 + offset, endpoint, news_endpoint)
                for derivative in derivatives:
                    _add(pieces[derivative], 6, tail[derivative])
            tensors[kind, rank] = pieces
            for derivative in derivatives:
                for grade, tensor in pieces[derivative].items():
                    for (ell, m), mode in project_stf(binary, rank, tensor, kind).items():
                        _add(result[derivative], (grade + offset, ell, m), mode)
    leading = leading_mode_coefficients(binary, position, velocity, ell_max=maximum + 2)
    for derivative in derivatives:
        for key, mode in leading[derivative].items():
            supplied_ranks = 7 if maximum == 7 and key[0] == key[1] - 2 else 6
            if key[0] <= maximum and key[1] > supplied_ranks:
                _add(result[derivative], key, mode)
    return result[0], result[1], tensors, current


def _project_endpoint(binary, kind, rank, grade, tensors, output):
    for power, tensor in tensors.items():
        for (ell, m), value in project_stf(binary, rank, tensor, kind).items():
            coefficients = output.setdefault((grade, ell, m), np.zeros(3, dtype=complex))
            coefficients[power] += value


def build_retained_memory(binary,canonical_history,radiative_scale,*,quadrupole_maximum_grade=8,
                          log_history=None,max_time=None):
    """Integrate the extra formal grades of the full unaveraged news product."""
    from .memory import memory_rate_polynomials, retained_memory_excess, RetainedMemory, _memory_weights

    moments = canonical_history.moments
    accepted = moments["orbital"][0].panels
    if max_time is not None:
        ends = np.array([panel.end for panel in accepted])
        last = np.searchsorted(ends,max_time,side="left")
        accepted = accepted[:last+1]
    origin, scale = accepted[0].start, 2 * radiative_scale
    endpoint = {}

    def sample(index, fraction, save_endpoint=False):
        current = {
            name: {
                grade: np.polynomial.polynomial.polyval(fraction, past.panels[index].coefficients)
                for grade, past in grades.items()
            }
            for name, grades in moments.items()
        }
        panel = accepted[index]
        time = panel.start + fraction * (panel.end - panel.start)
        position, velocity = current.pop("orbital")[0]
        return finite_radiation_available(
            binary, canonical_history, time, position, velocity, radiative_scale,
            news_only=True, sampled_jets=current, news_endpoint=endpoint if save_endpoint else False,
            log_history=log_history if log_history is not None and log_history.active(time) else None,
        )[1]

    snapshots = []
    total = _incoming_memory_excess(binary, canonical_history, radiative_scale, quadrupole_maximum_grade)
    first = sample(0, 0, save_endpoint=True)
    for index, panel in enumerate(accepted):
        middle, last = sample(index, 0.5), sample(index, 1)
        polynomial = {}
        for key, value in first.items():
            if key in endpoint:
                coefficients = np.zeros((3, 3), dtype=complex)
                coefficients[:, 0] = _interpolate((value, middle[key], last[key]))
                coefficients[0] += endpoint[key]
            else:
                coefficients = _interpolate((value, middle[key], last[key]))[:, None]
            polynomial[key] = coefficients
        rates, _ = retained_memory_excess(memory_rate_polynomials(binary, polynomial), {},
                                          quadrupole_maximum_grade=quadrupole_maximum_grade)
        snapshots.append((panel.start, panel.end, rates, total))
        weights = _memory_weights(panel.start, panel.end, origin, scale, 1)
        total = {key: total.get(key, 0) + np.sum(rate * weights[:, :rate.shape[1]])
                 for key, rate in rates.items()}
        first = last
    return RetainedMemory(tuple(snapshots), origin, scale)


def _incoming_memory_excess(binary, history, radiative_scale, quadrupole_maximum_grade):
    """Seed the higher formal memory with the same physical source approximation."""
    from .memory import null_memory_news, retained_memory_excess

    incoming = history.incoming
    if incoming is None:
        return {}
    accumulated = {}
    for time, weight in zip(*incoming.grid):
        local = replace(incoming, orbit=replace(incoming.orbit, cutoff=time))
        nested = replace(history, incoming=local)
        position, velocity = local.orbit.state_at(binary, time)
        canonical = local.join
        _, news, _, _ = finite_radiation_available(
            binary, nested, time, position, velocity, radiative_scale, news_only=True,
            sampled_jets=canonical,
        )
        rates, _ = retained_memory_excess(null_memory_news(binary, news), {},
                                          quadrupole_maximum_grade=quadrupole_maximum_grade)
        for key, value in rates.items():
            _add(accumulated, key, weight*value)
    return accumulated


def finite_quadrupole_grade_eight(
    binary, history, time, position, velocity, harmonic_scales, radiative_scale, *,
    hadamard_jet=None, leading_jet=None, prefix_cache=None, dynamical_tail_jet=None, log_history=None,
):
    """Compute the grade-eight quadrupole ingredients of a finite-history diagnostic.

    An independently obtained local 4PN Hadamard jet can be supplied when
    checking the full tensor. With no such jet the returned budget is partial.
    The current orbital state and the panels describe the same source time.
    """
    if history.maximum_absolute_grade != 7 or set(history.mass_grades) != {0, 2, 4}:
        raise ValueError("Provide the grade-seven canonical history with conserved mass grades 0, 2 and 4.")
    if harmonic_scales[0] != history.harmonic_scale:
        raise ValueError("The quadrupole and canonical histories require the same SH scale.")
    canonical = history.covered_moments()["M2"]
    leading_past = history.quadrupole_for_cubic_terms
    mass = history.mass_grades
    infrared_scale = history.infrared_scale
    current = leading_jet
    if current is None:
        current = canonical_waveform_jets(
            binary, position, velocity, history.harmonic_scale, infrared_scale, history.dipole_grades, 7,
        )["M2"][0]
    local = known_local_quadrupole_grade_eight(
        binary, position, velocity, harmonic_scales, infrared_scale,
    )[2:4]
    if hadamard_jet is not None:
        hadamard_jet = np.asarray(hadamard_jet, dtype=float)
        if hadamard_jet.shape != (2, 3, 3) or not np.isfinite(hadamard_jet).all():
            raise ValueError("Supply the second and third derivatives of the grade-eight Hadamard tensor.")
        local = local + hadamard_jet
    hereditary_source = nonlocal_quadrupole_grade_eight(
        binary, leading_past, time, position, velocity, infrared_scale, current,
        dynamical_tail_jet=dynamical_tail_jet,
        split_jet=None if log_history is None else log_history.leading,
    )[2:4]
    # The 2.5PN source-to-canonical conversion is included inside this ordinary tail.
    def split(grade):
        return None if log_history is None else partial(log_history.jet,"M2",grade,parts=2)

    ordinary = ordinary_tail(binary, canonical[5], time, 2, "mass", mass[0], radiative_scale,
                             split_jet=split(5))
    tailsquared = tail_of_tail(
        binary, canonical[2], time, 2, "mass", mass[0], radiative_scale, infrared_scale, split_jet=split(2),
    )
    tailsquared += 2 * mass[2] / mass[0] * tail_of_tail(
        binary, canonical[0], time, 2, "mass", mass[0], radiative_scale, infrared_scale, split_jet=split(0),
    )
    cubic = mpm_quadrupole_interactions(
        binary, leading_past, time, mass[0], history.dipole_grades[0],
        infrared_scale, radiative_scale, current=current, prefix_cache=prefix_cache, log_history=log_history,
    )[8]
    return {
        "canonical_local": local, "canonical_nonlocal": hereditary_source,
        "ordinary_tail": ordinary, "tail_of_tail": tailsquared, "cubic": cubic,
    }


def finite_radiation_grade_eight(
    binary, history, time, position, velocity, harmonic_scales, radiative_scale, hadamard_jet,
    *, prefix_cache=None, dynamical_tail_jet=None, log_history=None,
):
    """Assemble the 3.5PN mode hierarchy and an explicitly provided 4PN completion."""
    if hadamard_jet is None:
        raise ValueError("A grade-eight radiation calculation needs the complete local 4PN source jet.")
    strain, news, tensors, canonical = finite_radiation_available(
        binary,history,time,position,velocity,radiative_scale,log_history=log_history,prefix_cache=prefix_cache,
    )
    parts = finite_quadrupole_grade_eight(
        binary, history, time, position, velocity, harmonic_scales, radiative_scale,
        hadamard_jet=hadamard_jet, leading_jet=canonical["M2"][0], prefix_cache=prefix_cache,
        dynamical_tail_jet=dynamical_tail_jet,log_history=log_history,
    )
    quadrupole = sum(parts.values())
    for modes, graded, field in zip((strain, news), tensors["mass", 2], quadrupole):
        graded[8] = field
        modes.update({(8, ell, m): value for (ell, m), value in project_stf(binary, 2, field).items()})
    return strain, news, tensors, canonical, parts
