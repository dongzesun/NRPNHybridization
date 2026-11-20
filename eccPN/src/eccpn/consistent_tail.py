"""Derivative-consistent interpolation and integration of a causal tail primitive."""

from dataclasses import dataclass
from functools import cache
from math import comb, factorial

import numpy as np
from numpy.polynomial.legendre import leggauss

from .compact_moments import CompactMomentHistory, pack_stf, unpack_stf
from .history import batched_log_weights, integrate_log_panel


@dataclass(frozen=True)
class TailPrimitive:
    """Local Hermite panels with the same primitive jet at a shared endpoint."""

    cuts: np.ndarray
    coefficients: np.ndarray
    incoming: object

    @classmethod
    def from_canonical(cls, history, name, primitive_order, grade_weights, *, news_derivatives=1):
        moments = history.moments[name]
        panels = moments[0].panels
        cuts = np.array([panels[0].start, *(panel.end for panel in panels)])
        conditions = news_derivatives + 2
        orders = slice(primitive_order, primitive_order + conditions)
        first = 0
        last = 0
        for grade, weight in grade_weights.items():
            source = moments[grade].panels
            coefficients = np.stack([panel.coefficients[:, orders] for panel in source])
            first = first + weight * coefficients[:, 0]
            last = last + weight * coefficients.sum(axis=1)
        width = np.diff(cuts).reshape((-1,) + (1,) * (first.ndim - 2))
        powers = np.arange(conditions)
        dimensions = (1, conditions) + (1,) * (first.ndim-2)
        width_powers = width[:, None]**powers.reshape(dimensions)
        left = width_powers * first / np.array([factorial(i) for i in powers]).reshape(dimensions)
        endpoint, inverse = _endpoint_system(conditions)
        residual = width_powers*last - np.einsum("ij,pj...->pi...",endpoint,left)
        right = np.einsum("ij,pj...->pi...",inverse,residual)
        coefficients = np.concatenate((left,right),axis=1)
        incoming = None
        if history.incoming is not None:
            source = history.incoming
            values = sum(weight * source.sampled[name][grade][:, primitive_order]
                         for grade, weight in grade_weights.items())
            incoming = (*source.grid, values)
        return cls(cuts, coefficients, incoming)

    def integrate(self, time, lag, terms, *, nodes=None):
        """Return the tail and its derivative for common logarithmic kernel terms.

        Terms are (coefficient, scale, log power). The split is at fixed delay,
        and must fall in the covered interpolant. The remote past is integrated
        against the twice differentiated kernel.
        """
        split = time - lag
        if not self.cuts[0] < split < time <= self.cuts[-1] or self.incoming is None:
            raise ValueError("The split tail needs incoming data and a split inside the accepted panels.")
        index = np.searchsorted(self.cuts, split, side="right") - 1
        start, stop = self.cuts[index:index + 2]
        width = stop - start
        fraction = (np.longdouble(split)-start)/width
        primitive = self.coefficients[index].astype(np.longdouble)
        count = self.coefficients.shape[1]//2 - 1
        kernel = _kernel(terms,np.longdouble(time)-split,count+1)
        boundary = [np.polynomial.polynomial.polyval(fraction,_derive(primitive,order,width))
                    for order in range(count+1)]
        result = np.stack([sum(boundary[order+1-k]*kernel[k] for k in range(order+2))
                           for order in range(count)])

        times, weights, values = self.incoming
        result += np.tensordot(_kernel(terms,time-times,count+1)[2:] * weights, values, axes=1)
        result += self._far_integral(index,split,time,terms,nodes,count)
        end_index = np.searchsorted(self.cuts, time, side="left")
        for panel in range(index, end_index):
            a, b = self.cuts[panel:panel + 2]
            for output, order in enumerate(range(2,count+2)):
                coeff = _derive(self.coefficients[panel], order, b-a)
                lower = max(a, split)
                if lower != a:
                    coeff = _shift(coeff, (lower-a)/(b-a))
                for weight, scale, power in terms:
                    result[output] += weight * integrate_log_panel(lower, b, coeff, time, scale, power)
        return result

    def _far_integral(self,index,split,time,terms,nodes,count,*,first_order=2):
        starts = self.cuts[:index + 1]
        ends = np.minimum(self.cuts[1:index + 2], split)
        widths = np.diff(self.cuts[:index + 2])
        if nodes is not None:
            gauss_nodes,weights = leggauss(nodes)
        delay_start = np.longdouble(time)-starts
        delay_end = np.longdouble(time)-ends
        # Split broad panels geometrically in delay. This bounds the lag ratio
        # on each Gauss panel independently of the physical orbit's duration.
        sections = np.maximum(1, np.ceil(np.log2(delay_start/delay_end)).astype(int))
        total = np.zeros((count,) + self.coefficients.shape[2:],dtype=np.longdouble)
        for section in range(sections.max()):
            active = (sections > section) & (ends > starts)
            a = delay_start[active] * np.exp2(-section)
            b = np.maximum(a/2, delay_end[active])
            degree = self.coefficients.shape[1]-1
            if nodes is None:
                fractions = (delay_start[active]-a)/widths[active]
                span = (a-b)/widths[active]
                local = _inverse_log_panel_weights(a,b,terms,degree,first_order,count)
                panel_weights = np.zeros_like(local)
                for n in range(degree+1):
                    for j in range(n+1):
                        panel_weights[:,:,n] += comb(n,j)*fractions**(n-j)*span**j*local[:,:,j]
            else:
                half = (a-b)/2
                delay = b[:,None]+half[:,None]*(1+gauss_nodes)
                fraction = (delay_start[active,None]-delay)/widths[active,None]
                kernels = _kernel(terms,delay,first_order+count-1)[first_order:]*(half[:,None]*weights)
                powers = fraction[:,:,None]**np.arange(degree+1)
                panel_weights = np.einsum("kij,ijn->kin",kernels,powers)
            total += np.tensordot(panel_weights, self.coefficients[:index + 1][active], axes=2)
        return total


def _derive(coefficients, order, width):
    result = coefficients
    for _ in range(order):
        shape = (-1,) + (1,)*(result.ndim-1)
        result = result[1:] * np.arange(1, len(result)).reshape(shape) / width
    return result


def _shift(coefficients, fraction):
    remainder = 1-fraction
    return np.stack([
        remainder**degree * sum(comb(old, degree) * fraction**(old-degree) * coefficients[old]
                               for old in range(degree, len(coefficients)))
        for degree in range(len(coefficients))
    ])


@cache
def _endpoint_system(conditions):
    derivatives = np.array([[factorial(degree)//factorial(degree-order) if degree>=order else 0
                             for degree in range(2*conditions)] for order in range(conditions)])
    return derivatives[:,:conditions],np.linalg.inv(derivatives[:,conditions:])


@cache
def _inverse_series(order, degree):
    # (1-z)^(-order), and its derivative with respect to order.
    coefficients = []
    derivatives = []
    logarithmic_slope = np.longdouble(0)
    for n in range(112):
        binomial = np.longdouble(comb(n+order-1,n))/(n+degree+1)
        coefficients.append(binomial)
        derivatives.append(binomial*logarithmic_slope)
        logarithmic_slope += np.longdouble(1)/(order+n)
    return np.array(coefficients),np.array(derivatives)


def _divided_series(values,ratio):
    # Bound each truncated positive series by its geometric tail.
    # The extra polynomial factors are at most degree six in this use.
    result = np.empty_like(ratio)
    small = ratio<=1/64
    medium = (ratio>1/64)&(ratio<=1/8)
    for mask,terms in ((small,16),(medium,40),(~(small|medium),len(values))):
        result[mask] = np.polynomial.polynomial.polyval(ratio[mask],values[:terms])
    return result


def _inverse_log_panel_weights(start,stop,terms,degree,first,count):
    width = start-stop
    ratio = width/start
    weights = np.zeros((count,len(start),degree+1),dtype=np.longdouble)
    for output,order in enumerate(range(first,first+count)):
        harmonic = sum(np.longdouble(1)/i for i in range(1,order))
        kernel_factor = 0
        log_factor = 0
        for coefficient,scale,power in terms:
            if power == 1:
                kernel_factor += coefficient
            elif power == 2:
                kernel_factor += 2*coefficient*(np.log(start/scale)-harmonic)
                log_factor -= 2*coefficient
        factor = (-1)**(order-1)*factorial(order-1)*width/start**order
        for j in range(degree+1):
            series,derivative = _inverse_series(order,j)
            value = kernel_factor*_divided_series(series,ratio)
            if log_factor:
                value += log_factor*_divided_series(derivative,ratio)
            weights[output,:,j] = factor*value
    return weights


def _kernel(terms,delay,highest):
    value = np.zeros((highest+1,) + np.shape(delay),dtype=np.result_type(delay,float))
    for weight,scale,power in terms:
        log = np.log(delay/scale)
        value[0] += weight*log**power
        if power == 0:
            continue
        harmonic = np.array(0,dtype=value.dtype)
        for order in range(1,highest+1):
            factor = weight*(-1)**(order-1)*factorial(order-1)/delay**order
            value[order] += factor if power==1 else 2*factor*(log-harmonic)
            harmonic += np.asarray(1,dtype=value.dtype)/order
    return value


def _near_integral(cuts,coefficients,index,split,time,terms):
    end = np.searchsorted(cuts,time,side="left")
    first_start,first_end = cuts[index:index+2]
    clipped = _shift(coefficients[index],(split-first_start)/(first_end-first_start))
    result = np.zeros(coefficients.shape[2:])
    starts,ends = cuts[index+1:end],cuts[index+2:end+1]
    for coefficient,scale,power in terms:
        result += coefficient*integrate_log_panel(split,first_end,clipped,time,scale,power)
        if starts.size:
            weights = batched_log_weights(starts,ends,time,scale,power,coefficients.shape[1])
            result += coefficient*np.tensordot(weights,coefficients[index+1:end],axes=2)
    return result


@dataclass(frozen=True)
class SeparatedLogJet:
    """Exact near log weights and analytic integrations by parts in the far region.

    The same formal canonical jet supplies the split boundary and the near
    convolution. Higher time derivatives of the far kernel are taken before
    quadrature, avoiding cumulative offsets from high local derivatives.
    """

    cuts: np.ndarray
    jets: np.ndarray
    incoming_times: np.ndarray
    incoming_weights: np.ndarray
    incoming_values: np.ndarray
    rank: int | None = None

    @classmethod
    def from_canonical(cls, moment):
        compact = isinstance(moment.finite,CompactMomentHistory)
        rank = moment.finite.rank if compact else None
        panels = moment.finite.packed.panels if compact else moment.panels
        cuts = np.array([panels[0].start,*(panel.end for panel in panels)])
        times,weights = moment.past.grid
        values = moment.past.sampled[moment.name][moment.grade]
        if rank is not None:
            values = pack_stf(values,rank)
        return cls(cuts,np.stack([panel.coefficients for panel in panels]),times,weights,values,rank)

    def integrate(self,time,lag,terms,order,count,parts=4,*,nodes=None):
        """Integrate the jet orders [order, order + count) with a common split."""
        split = time-lag
        if not self.cuts[0]<split<time<=self.cuts[-1] or not 1<=parts<=order:
            raise ValueError("The far integration must split inside the committed history.")
        index = np.searchsorted(self.cuts,split,side="right")-1
        start,stop = self.cuts[index:index+2]
        fraction = (np.longdouble(split)-start)/(stop-start)
        boundary = np.polynomial.polynomial.polyval(fraction,self.jets[index].astype(np.longdouble))
        highest = parts+count-1
        kernel = _kernel(terms,np.longdouble(time)-split,highest)
        base = order-parts
        result = np.stack([sum(boundary[order+j-1-k].astype(np.longdouble)*kernel[k] for k in range(parts+j))
                           for j in range(count)])
        far = TailPrimitive(self.cuts,self.jets[:,:,base],None)
        result += far._far_integral(index,split,time,terms,nodes,count,first_order=parts)
        delay = time-self.incoming_times
        remote_weights = _kernel(terms,delay,highest)[parts:]*self.incoming_weights
        result += np.tensordot(remote_weights,self.incoming_values[:,base],axes=1)
        coefficients = self.jets[:,:,order:order+count]
        result += _near_integral(self.cuts,coefficients,index,split,time,terms)
        result = np.asarray(result,dtype=float)
        return unpack_stf(result,self.rank) if self.rank is not None else result


@dataclass(frozen=True)
class SeparatedLogRate:
    """Analytic rate of the same tabulated continuous product used in a log integral."""

    cuts: np.ndarray
    products: np.ndarray
    incoming_times: np.ndarray
    incoming_weights: np.ndarray
    incoming_values: np.ndarray

    def value(self,time,terms):
        count = np.searchsorted(self.cuts[:-1],time,side="left")
        starts,ends = self.cuts[:count],self.cuts[1:count+1]
        delay = time-self.incoming_times
        value = np.zeros(self.products.shape[2:])
        for coefficient,scale,power in terms:
            weights = batched_log_weights(starts,ends,time,scale,power,self.products.shape[1])
            value += coefficient*np.tensordot(weights,self.products[:count],axes=2)
            weights = self.incoming_weights*np.log(delay/scale)**power
            value += coefficient*np.tensordot(weights,self.incoming_values,axes=1)
        return value

    def rate(self,time,lag,terms,*,nodes=None):
        split = time-lag
        if not self.cuts[0]<split<time<=self.cuts[-1]:
            raise ValueError("The product rate needs a split inside the committed history.")
        index = np.searchsorted(self.cuts,split,side="right")-1
        start,stop = self.cuts[index:index+2]
        fraction = (np.longdouble(split)-start)/(stop-start)
        boundary = np.polynomial.polynomial.polyval(fraction,self.products[index].astype(np.longdouble))
        result = boundary.astype(np.longdouble)*_kernel(terms,np.longdouble(time)-split,0)[0]
        far = TailPrimitive(self.cuts,self.products,None)
        result += far._far_integral(index,split,time,terms,nodes,1,first_order=1)[0]
        delay = time-self.incoming_times
        weight = _kernel(terms,delay,1)[1]*self.incoming_weights
        result += np.tensordot(weight,self.incoming_values,axes=1)
        end = np.searchsorted(self.cuts,time,side="left")
        for panel in range(index,end):
            a,b = self.cuts[panel:panel+2]
            coeff = _derive(self.products[panel],1,b-a)
            lower = max(a,split)
            if lower!=a:
                coeff = _shift(coeff,(lower-a)/(b-a))
            for coefficient,scale,power in terms:
                result += coefficient*integrate_log_panel(lower,b,coeff,time,scale,power)
        return np.asarray(result,dtype=float)


class RadiativeLogHistory:
    """Per-radiation-call cache of source and product histories for the split backend."""

    def __init__(self, history, lag):
        if history.incoming is None:
            raise ValueError("The split radiative backend requires a matched incoming waveform.")
        self.moments = history.covered_moments()
        self.lag = lag
        self._linear = {}
        self._products = {}
        self._time = None
        self._integrals = {}

    def active(self, time):
        return time-self.lag > self.moments["M2"][0].start

    def jet(self,name,grade,time,terms,order,count,*,parts=4):
        key = (name,grade)
        if key not in self._linear:
            self._linear[key] = SeparatedLogJet.from_canonical(self.moments[name][grade])
        if time!=self._time:
            self._time = time
            self._integrals.clear()
        request = (*key,tuple(terms),order,count,parts)
        if request not in self._integrals:
            self._integrals[request] = self._linear[key].integrate(time,self.lag,terms,order,count,parts)
        return self._integrals[request].copy()

    def leading(self, time, terms, order, count):
        return self.jet("M2",0,time,terms,order,count)

    def product(self,time,left,right,terms):
        from .radiative_quadrupole import _contract, _convolve_coefficients

        key = (left,right)
        if key not in self._products:
            name,grade = "M2",0
            if (name,grade) not in self._linear:
                self._linear[name,grade] = SeparatedLogJet.from_canonical(self.moments[name][grade])
            source = self._linear[name,grade]
            products = np.stack([_convolve_coefficients(panel[:,left],panel[:,right])
                                 for panel in source.jets])
            incoming = _contract(source.incoming_values[:,left],source.incoming_values[:,right])
            self._products[key] = SeparatedLogRate(source.cuts,products,source.incoming_times,
                                                   source.incoming_weights,incoming)
        past = self._products[key]
        return np.stack((past.value(time,terms),past.rate(time,self.lag,terms)))
