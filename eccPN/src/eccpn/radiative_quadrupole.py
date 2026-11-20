"""Explicit 2.5PN, 3PN and cubic 4PN MPM mass quadrupole interactions."""

from bisect import bisect_right
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad_vec

from .compressed_history import CompressedPanelHistory
from .history import SourcePanel, integrate_log_panel


def _stf(tensor):
    return (tensor + np.swapaxes(tensor, -2, -1)) / 2 - (
        np.trace(tensor, axis1=-2, axis2=-1)[..., None, None] * np.eye(3) / 3
    )


def _contract(left, right):
    return _stf(left @ right)


def _spin_contract(spin, moment):
    return _stf(np.cross(spin[..., None, :], moment).swapaxes(-2, -1))


def _convolve_coefficients(first, second):
    products = np.einsum("aij,bjk->abik", first, second)
    degree = len(first) + len(second) - 1
    result = np.zeros((degree, 3, 3))
    for index, row in enumerate(products):
        result[index:index + len(second)] += row
    return _stf(result)


def _kernel_value(lag, terms):
    return sum(coefficient * np.log(lag / scale)**power for coefficient, scale, power in terms)


def _linear_integral(history, time, derivative, kernel, split_jet=None):
    if split_jet is not None:
        return split_jet(time,kernel,derivative,2)
    value, news = sum(
        (coefficient * history.integrate_jet(time, scale, power, slice(derivative, derivative+2))
         for coefficient, scale, power in kernel),
        np.zeros((2, 3, 3)),
    )
    news = news + history.start_jump(derivative) * _kernel_value(time - history.start, kernel)
    return value, news


def _product_coefficients(panel, left, right):
    moments = panel.coefficients
    pair = _convolve_coefficients(moments[:, left], moments[:, right])
    derivative = (_convolve_coefficients(moments[:, left + 1], moments[:, right])
                  + _convolve_coefficients(moments[:, left], moments[:, right + 1]))
    return np.stack((pair, derivative), axis=1)


def _quadratic_integral(history, time, left, right, kernel, products=None):
    result = np.zeros((2,3,3))
    for index, panel in enumerate(history.panels):
        if panel.start >= time:
            break
        if products is None:
            polynomial = _product_coefficients(panel, left, right)
        else:
            cache = products.setdefault((left, right), [])
            if len(cache) <= index:
                cache.append(_product_coefficients(panel, left, right))
            polynomial = cache[index]
        for coefficient,scale,power in kernel:
            result += coefficient*integrate_log_panel(panel.start,panel.end,polynomial,time,scale,power)
    initial = history.value(history.start)
    jump = _contract(initial[left],initial[right])
    remote = getattr(history,"remote",None)
    if remote is not None:
        for coefficient,scale,power in kernel:
            result[0] += coefficient*remote.product_integral(time,scale,power,left,right)
            result[1] += coefficient*(remote.product_integral(time,scale,power,left+1,right)
                                       +remote.product_integral(time,scale,power,left,right+1))
        jump -= _contract(remote.value(history.start,left),remote.value(history.start,right))
    result[1] += jump*_kernel_value(time-history.start,kernel)
    return result


def _instantaneous_product(moments, left, right):
    return np.stack((
        _contract(moments[left], moments[right]),
        _contract(moments[left + 1], moments[right]) + _contract(moments[left], moments[right + 1]),
    ))


@dataclass(frozen=True)
class _NestedSource:
    history: object
    scale: float
    offset: float
    rtol: float
    atol: float
    inner: object

    @classmethod
    def prepare(cls, history, scale, offset, rtol, atol):
        panels = history.panels
        inner = None
        if len(panels) >= 32 and all(len(panel.coefficients) <= 3 for panel in panels):
            projected = (SourcePanel(panel.start, panel.end, panel.coefficients[:, 4]) for panel in panels)
            inner = CompressedPanelHistory(history.start, (3,3), "ordered-tail-inner", None).append(*projected)
        return cls(history, scale, offset, rtol, atol, inner)

    def rate(self, time):
        history, scale, offset = self.history, self.scale, self.offset
        source = history.value(time)[4]
        if self.inner is None:
            lagged = history.integrate_jet(time, scale, 1, 4)
            lagged += offset * history.integrate_jet(time, scale, 0, 4)
        else:
            lagged = self.inner.integral(time, scale*np.exp(-offset))
            remote = getattr(history, "remote", None)
            if remote is not None:
                lagged += remote.linear_integral(time, scale, 1, 4)
                lagged += offset*remote.linear_integral(time, scale, 0, 4)
        return _contract(source, lagged)

    def panel_integral(self, index, end):
        contribution, _, info = quad_vec(
            self.rate, self.history.panels[index].start, end,
            epsabs=self.atol, epsrel=self.rtol, full_output=True,
        )
        if not info.success:
            raise RuntimeError(f"The ordered tail-of-memory integral failed: {info.message}")
        return contribution


def _nested_full_panels(source, indices):
    return [source.panel_integral(index, source.history.panels[index].end) for index in indices]


def _new_nested_prefix(history, scale, offset, rtol, atol):
    remote = getattr(history, "remote", None)
    initial = np.zeros((3,3)) if remote is None else remote.nested_remote(scale, offset, rtol=rtol, atol=atol)
    return [initial], tuple(panel.end for panel in history.panels), _NestedSource.prepare(history, scale, offset, rtol, atol)


def _extend_nested_prefix(prefix, time, executor=None):
    from .parallel import ordered_chunks

    cumulative, ends, source = prefix
    completed = bisect_right(ends, time)
    indices = range(len(cumulative) - 1, completed)
    for rows in ordered_chunks(executor, source, _nested_full_panels, indices):
        for contribution in rows:
            cumulative.append(cumulative[-1] + contribution)


def _nested_memory(history, time, scale, rtol, atol, offset=-1613 / 270, prefix_cache=None):
    key = (id(history), scale, offset, rtol, atol)
    if prefix_cache is None:
        prefix_cache = {}
    if key not in prefix_cache:
        prefix_cache[key] = _new_nested_prefix(history, scale, offset, rtol, atol)
    prefix = prefix_cache[key]
    _extend_nested_prefix(prefix, time)
    cumulative, ends, source = prefix
    completed = bisect_right(ends, time)
    nested = cumulative[completed]
    if completed != len(ends) and history.panels[completed].start < time:
        nested = nested + source.panel_integral(completed, time)
    return np.stack((nested, source.rate(time)))


def mpm_quadrupole_interactions(
    binary, history, time, canonical_mass, current_dipole, infrared_scale, radiative_scale,
    *, rtol=1e-9, atol=1e-18, current=None, prefix_cache=None, log_history=None, cubic_only=False,
):
    """Return U_ij and its news by explicit grade for a finite jet history.

    Each panel represents derivatives zero through eight of one canonical
    mass quadrupole. The input current dipole is the conserved orbital
    angular momentum for nonspinning bodies. The radiative scale is a time.
    """
    if (history.source_shape not in ((9, 3, 3), (10, 3, 3)) or not history.start < time <= history.end):
        raise ValueError("This quadrupole interaction requires a covered finite history of at least nine STF jets.")
    if (not np.isfinite([canonical_mass, infrared_scale, radiative_scale, rtol, atol]).all()
            or min(canonical_mass, infrared_scale, radiative_scale, rtol, atol) <= 0):
        raise ValueError("Canonical mass, scales and tolerances must be finite and positive.")
    spin = np.asarray(current_dipole, dtype=float)
    if spin.shape != (3,) or not np.isfinite(spin).all():
        raise ValueError("The conserved current dipole must be a finite Cartesian vector.")
    current = history.value(time) if current is None else np.asarray(current, dtype=float)
    if current.shape != history.source_shape or not np.isfinite(current).all():
        raise ValueError("Supply all Cartesian quadrupole derivatives at the current cut.")
    scale_ir = 2 * infrared_scale / binary.c
    scale_rad = 2 * radiative_scale
    constant = ((1, scale_ir, 0),)
    split_jet = None if log_history is None else log_history.leading

    lower = {}
    if not cubic_only:
        # QC2 Eq. (2.5), the quadratic 2.5PN mass quadrupole.
        quadratic = (
            -(2 / 7) * _quadratic_integral(history, time, 3, 3, constant)
            + (1 / 7) * _instantaneous_product(current, 5, 0)
            - (5 / 7) * _instantaneous_product(current, 4, 1)
            - (2 / 7) * _instantaneous_product(current, 3, 2)
            - (1 / 3) * _spin_contract(spin, current[4:6])
        )
        tail_tail = np.stack(_linear_integral(history, time, 5, (
            (1, scale_rad, 2), (11 / 6, scale_rad, 1),
            (-107 / 105, scale_ir, 1), (124627 / 44100, scale_ir, 0),
        )))
        lower[5] = binary.G / binary.c**5 * quadratic
        lower[6] = 2 * (binary.G**2 * canonical_mass) * canonical_mass / binary.c**6 * tail_tail

    # QC2 Eq. (2.7), the generic 4PN cubic interactions.
    cubic = _nested_memory(history, time, scale_ir, rtol, atol, prefix_cache=prefix_cache)
    products = (None if prefix_cache is None else
                prefix_cache.setdefault(("quadratic_products", id(history)), {}))
    for left, right, coefficient, rad_log in (
        (3, 4, -5 / 2, 3 / 2), (2, 5, -3, 11 / 12),
        (1, 6, -5 / 2, 3 / 10), (0, 7, -1, -1 / 4),
    ):
        kernel = ((1,scale_ir,1),(rad_log,scale_rad,1))
        if log_history is None:
            integral = _quadratic_integral(history, time, left, right, kernel, products)
        else:
            integral = log_history.product(time,left,right,kernel)
        cubic += coefficient*integral
    for left, right, coefficient, offset in (
        (2, 5, -2, 27521 / 5040), (1, 6, -5 / 2, 15511 / 3150), (0, 7, 1 / 2, -6113 / 756),
    ):
        past, past_news = _linear_integral(history, time, right, ((1,scale_ir,1),(offset,scale_ir,0)),split_jet)
        cubic += coefficient * np.stack((
            _contract(current[left], past),
            _contract(current[left + 1], past) + _contract(current[left], past_news),
        ))
    current_tail = _spin_contract(spin, np.stack(_linear_integral(history, time, 6, (
        (1, scale_rad, 1), (2, scale_ir, 1), (1223 / 1890, scale_ir, 0),
    ),split_jet)))
    gm_factor = binary.G**2 * canonical_mass
    lower[8] = gm_factor / binary.c**8 * ((8 / 7) * cubic - (2 / 3) * current_tail)
    return lower
