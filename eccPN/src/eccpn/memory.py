"""Finite graded Bondi-news products for the null memory source."""

from dataclasses import dataclass
from functools import cache, cached_property
from importlib.resources import files
from math import comb, sqrt

import numpy as np
from scipy.sparse import csr_matrix

from .history import integrate_log_panel
from .projection import project_stf, stf_tensor


@cache
def _couplings():
    path = files("eccpn").joinpath("generated", "memory_couplings.npz")
    with path.open("rb") as source, np.load(source, allow_pickle=False) as arrays:
        indices, values = arrays["indices"], arrays["values"]
    result = {}
    for (ell1, m1, ell2, m2, ell, m), value in zip(indices, values):
        key = (int(ell1), int(m1), int(ell2), int(m2))
        result.setdefault(key, []).append((int(ell), int(m), value))
    return result


def angular_news_products(news, max_grade=7):
    """Project products of H-news to scalar harmonics by explicit 1/c grade.

    Input dictionaries use keys (grade, ell, m). Signed magnetic indices are
    required. The result has the same key layout, with grades added before
    truncation. This product is proportional to the instantaneous Bondi flux.
    """
    return _angular_products(news, max_grade, lambda first, second: first * np.conj(second))


def _angular_products(news, max_grade, multiply):
    if max_grade not in range(8):
        raise ValueError("The generated angular table covers product grades 0 to 7.")
    grouped = {}
    for (grade, ell, m), series in news.items():
        if ell not in range(2, 10) or abs(m) > ell or grade not in range(ell - 2, 8):
            raise ValueError("Input news mode is outside the supported fixed-history grade inventory.")
        grouped.setdefault((ell, m), []).append((grade, np.asarray(series)))
    products = {}
    for (ell1, m1), first in grouped.items():
        for (ell2, m2), second in grouped.items():
            couplings = _couplings().get((ell1, m1, ell2, m2), ())
            if not couplings:
                continue
            for grade1, mode1 in first:
                for grade2, mode2 in second:
                    grade = grade1 + grade2
                    if grade > max_grade:
                        continue
                    amplitude = multiply(mode1, mode2)
                    for ell, m, weight in couplings:
                        key = (grade, ell, m)
                        products[key] = products.get(key, 0) + weight * amplitude
    return products


def null_memory_news(binary, news, max_grade=7):
    """Compute the raw integrated-memory source dH_mem/du on each native cut."""
    result = {}
    for (grade, ell, m), product in angular_news_products(news, max_grade).items():
        spin_factor = 1 / sqrt((ell - 1) * ell * (ell + 1) * (ell + 2))
        result[grade + 5, ell, m] = (spin_factor / binary.c) * product
    return result


def integrate_null_memory_panel(binary, news, duration, fraction=1.0, max_grade=7):
    """Integrate the graded news flux source for polynomial H-news on one panel."""
    if not np.isfinite([duration, fraction]).all() or duration <= 0 or not 0 <= fraction <= 1:
        raise ValueError("A memory panel requires a positive width and fraction in [0, 1].")
    rates = memory_rate_polynomials(binary, news, max_grade)
    return {key: _integrate_rate(values, duration, fraction) for key, values in rates.items()}


@cache
def _polynomial_couplings(keys, max_grade):
    """Compile the same ordered angular products for repeated panels."""
    if max_grade not in range(8):
        raise ValueError("The generated angular table covers product grades 0 to 7.")
    pairs, rows, columns, weights, targets = [], [], [], [], {}
    for first, (grade_a, ell_a, m_a) in enumerate(keys):
        if ell_a not in range(2, 10) or abs(m_a) > ell_a or grade_a not in range(ell_a - 2, 8):
            raise ValueError("Input news mode is outside the supported fixed-history grade inventory.")
        for second, (grade_b, ell_b, m_b) in enumerate(keys):
            grade = grade_a + grade_b
            if grade > max_grade:
                continue
            terms = _couplings().get((ell_a, m_a, ell_b, m_b))
            if terms:
                for ell, m, weight in terms:
                    key = (grade + 5, ell, m)
                    row = targets.setdefault(key, len(targets))
                    rows.append(row)
                    columns.append(len(pairs))
                    weights.append(weight / sqrt((ell - 1) * ell * (ell + 1) * (ell + 2)))
                pairs.append((first, second))
    pair_indices = np.asarray(pairs, dtype=np.intp).reshape((-1, 2)).T
    matrix = csr_matrix((weights, (rows, columns)), shape=(len(targets), len(pairs)))
    return tuple(targets), pair_indices, matrix


def memory_rate_polynomials(binary, news, max_grade=7):
    """Return the analytic rate polynomials on the common unit panel."""
    keys = tuple(sorted(news))
    targets, (first, second), matrix = _polynomial_couplings(keys, max_grade)
    if not targets:
        return {}
    fields = [np.asarray(news[key]) for key in keys]
    has_logs = any(field.ndim == 2 for field in fields)
    fields = [field[:, None] if field.ndim == 1 else field for field in fields]
    n_times = max(len(field) for field in fields)
    n_logs = max(field.shape[1] for field in fields)
    polynomial = np.zeros((len(keys), n_times, n_logs), dtype=np.result_type(*fields))
    for index, field in enumerate(fields):
        polynomial[index, :field.shape[0], :field.shape[1]] = field
    left, right = polynomial[first], polynomial[second].conj()
    products = np.zeros((len(first), 2 * n_times - 1, 2 * n_logs - 1), dtype=polynomial.dtype)
    for t in range(n_times):
        for log in range(n_logs):
            coefficient = left[:, t, log]
            if np.any(coefficient):
                products[:, t:t + n_times, log:log + n_logs] += coefficient[:, None, None] * right
    if np.any(products[:, :, 3:]):
        raise ValueError("The retained news product exceeds the squared-log integration budget.")
    if has_logs:
        products = np.pad(products[:, :, :3], ((0, 0), (0, 0), (0, max(0, 3 - products.shape[2]))))
    else:
        products = products[:, :, 0]
    rates = matrix @ products.reshape((len(first), -1)) / binary.c
    return dict(zip(targets, rates.reshape((len(targets), *products.shape[1:]))))


def _integrate_rate(coefficients, duration, fraction):
    orders = np.arange(1, len(coefficients) + 1)
    return duration * np.sum(coefficients * fraction**orders / orders)


def _memory_weights(start, stop, origin, scale, fraction):
    """Weights for up to fourth degree and squared logs from the fixed initial cut."""
    if fraction == 0:
        return np.zeros((5, 3))
    cut = start + fraction * (stop - start)
    unit = np.eye(5)
    reversed_weights = np.column_stack([
        integrate_log_panel(-cut, -start, unit, -origin, scale, power) for power in range(3)
    ])
    return np.array([
        fraction**order * sum((-1)**j * comb(order, j) * reversed_weights[j] for j in range(order + 1))
        for order in range(5)
    ])


@dataclass(frozen=True)
class RetainedMemory:
    """Polynomial excess memory on shared finite panels, with cumulative strains."""

    panels: tuple
    origin: float
    scale: float

    @cached_property
    def _ends(self):
        return np.array([panel[1] for panel in self.panels])

    def value(self, time):
        if not self.panels[0][0] <= time <= self.panels[-1][1]:
            raise ValueError("The represented memory must cover the requested time.")
        if time == self.origin:
            raise ValueError("The finite-start memory rate requires a cut strictly after the initial cut.")
        start,stop,rates,prefix = self.panels[np.searchsorted(self._ends,time,side="left")]
        fraction = (time - start) / (stop - start)
        weights = _memory_weights(start, stop, self.origin, self.scale, fraction)
        lag_log = np.log((time - self.origin) / self.scale)
        strain = {key: prefix.get(key, 0) + np.sum(rate * weights)
                  for key, rate in rates.items()}
        news = {key: np.polynomial.polynomial.polyval(lag_log, np.polynomial.polynomial.polyval(fraction, rate))
                for key, rate in rates.items()}
        return strain, news


def _boundary_tensor(first, second, rank):
    if rank == 2:
        return stf_tensor(np.einsum("...ab,...ijab->...ij", first, second), rank)
    if rank == 4:
        return stf_tensor(np.einsum("...ai,...jkla->...ijkl", first, second), rank)
    return stf_tensor(np.einsum("...ij,...klmn->...ijklmn", first, second), rank)


def memory_boundary(binary, canonical_jets):
    """Return the formal grade-seven raw-to-MPM boundary strain and its news."""
    first = canonical_jets["M2"][0]
    second = canonical_jets["M4"][0]
    strain, news = {}, {}
    for rank, coefficient, power in ((2, 5 / 756, 7), (4, -12 / 55, 5), (6, 15 / 14, 3)):
        scale = binary.G * coefficient / binary.c**power
        tensor = scale * _boundary_tensor(first[3], second[4], rank)
        derivative = scale * (
            _boundary_tensor(first[4], second[4], rank) + _boundary_tensor(first[3], second[5], rank)
        )
        for (ell, m), value in project_stf(binary, rank, tensor).items():
            strain[7, ell, m] = value
        for (ell, m), value in project_stf(binary, rank, derivative).items():
            news[7, ell, m] = value
    return strain, news


def raw_quadrupole_overlap_grade_eight(binary, history, time, mass, radiative_scale, *, rtol=1e-9, atol=1e-18):
    """Rank-two raw memory overlap from Newtonian news times its ordinary tail.

    This formula uses integration by parts in the outer radiation-time
    integral. It preserves the derivative of the finite-start tail kernel
    and assumes the supplied jets are derivatives of a common source.
    """
    from .radiation_kernels import ordinary_tail
    from .radiative_quadrupole import _contract, _nested_memory

    if (history.source_shape not in ((9, 3, 3), (10, 3, 3))
            or not history.start < time <= history.end):
        raise ValueError("Supply a covered finite quadrupole history through the required jet orders.")
    if not np.isfinite([mass, radiative_scale, rtol, atol]).all() or min(mass, radiative_scale, rtol, atol) <= 0:
        raise ValueError("Mass, kernel scale and tolerances must be finite and positive.")
    quadrupole_news = history.value(time)[3]
    tail, tail_news = ordinary_tail(binary, history, time, 2, "mass", mass, radiative_scale)
    integrated, _ = _nested_memory(history, time, 2 * radiative_scale, rtol, atol, offset=11 / 12)
    endpoint = _quadrupole_tail_join(binary, history, mass, radiative_scale)
    return np.stack((
        8 * binary.G**2 * mass / (7 * binary.c**8) * integrated
        - 4 * binary.G / (7 * binary.c**5) * (_contract(quadrupole_news, tail) - endpoint),
        -4 * binary.G / (7 * binary.c**5) * _contract(quadrupole_news, tail_news),
    ))


def quadrupole_tail_memory_boundary(binary, history, time, mass, infrared_scale, radiative_scale):
    """Convert the QC2 nested term to the rank-two raw news-product memory.

    This is the boundary for the nested term of QC2 Eq. (2.7). Other
    terms of that equation remain in the fixed-history radiative tensor.
    """
    from .radiation_kernels import ordinary_tail
    from .radiative_quadrupole import _contract

    tail, tail_news = ordinary_tail(binary, history, time, 2, "mass", mass, radiative_scale)
    third, fourth = history.value(time)[3:5]
    difference = third - history.start_jump(3)
    lag_constant = np.log(infrared_scale / (binary.c * radiative_scale)) + 11 / 12 + 1613 / 270
    cut = 8 * binary.G**2 * mass / (7 * binary.c**8) * lag_constant
    product = 4 * binary.G / (7 * binary.c**5)
    return np.stack((
        cut / 2 * _contract(difference, difference)
        - product * (_contract(third, tail) - _quadrupole_tail_join(binary, history, mass, radiative_scale)),
        cut * _contract(difference, fourth)
        - product * (_contract(fourth, tail) + _contract(third, tail_news)),
    ))


def _quadrupole_tail_join(binary, history, mass, radiative_scale):
    from .radiation_kernels import ordinary_tail
    from .radiative_quadrupole import _contract

    if getattr(history, "remote", None) is None:
        return 0.0
    tail, _ = ordinary_tail(binary, history, history.start, 2, "mass", mass, radiative_scale, endpoint=False)
    return _contract(history.start_jump(3), tail)


def retained_memory_excess(raw_strain, raw_news, *, quadrupole_maximum_grade=8, other_maximum_grade=7):
    """Extract raw terms above each mode's fixed-history MPM budget."""
    if quadrupole_maximum_grade not in (7, 8) or other_maximum_grade not in range(5, 8):
        raise ValueError("Choose the absolute grade of the higher-mode MPM memory.")
    return tuple(
        {key: value for key, value in source.items()
         if key[0] > (quadrupole_maximum_grade if key[1] == 2 else other_maximum_grade)}
        for source in (raw_strain, raw_news)
    )
