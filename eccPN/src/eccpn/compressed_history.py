"""Immutable accepted-panel tree for the causal log and log-squared kernels."""

from dataclasses import dataclass
from functools import cached_property
from math import comb

import numpy as np


_ORDER = 24
_RATIO = 1/4
_DEGREE = np.arange(_ORDER+1)
_BASE_WEIGHTS = np.array([
    [sum(comb(p,k)*2/(n+k+1) for k in range(p+1) if (n+k)%2 == 0)/2**(p+1)
     for p in range(3)] for n in range(_ORDER+1)
])
_MERGE_COEFFICIENTS = np.array([
    [comb(n,k) if k<=n else 0 for k in _DEGREE] for n in _DEGREE
],dtype=float)
_SHIFT_ORDERS = np.maximum(0,_DEGREE[:,None]-_DEGREE)
_HARMONIC = np.r_[0.,0.,np.cumsum(1/np.arange(1,_ORDER))]
_INVERSE_WEIGHTS = {power: np.array([comb(n+power-1,power-1) for n in _DEGREE])
                    for power in (4,5)}


@dataclass(frozen=True)
class _Block:
    start: float
    end: float
    count: int
    moments: np.ndarray
    members: tuple

    @classmethod
    def leaf(cls,panel):
        half_width = (panel.end-panel.start)/2
        data = panel.coefficients
        moments = 2*half_width*np.tensordot(_BASE_WEIGHTS[:,:len(data)],data,axes=1)
        moments.flags.writeable = False
        return cls(panel.start,panel.end,1,moments,(panel,))

    @classmethod
    def merge(cls,left,right):
        # Transform each child's normalized coordinates to the parent.
        total_width = right.end-left.start
        left_fraction = (left.end-left.start)/total_width
        right_fraction = (right.end-right.start)/total_width
        weights = [
            _MERGE_COEFFICIENTS*(offset**_SHIFT_ORDERS)*(fraction**_DEGREE)[None,:]
            for offset,fraction in ((-right_fraction,left_fraction),(left_fraction,right_fraction))
        ]
        moments = np.tensordot(weights[0],left.moments,axes=1)+np.tensordot(weights[1],right.moments,axes=1)
        moments.flags.writeable = False
        return cls(left.start,right.end,left.count+right.count,moments,(left,right))

    def contracted(self,time,scale,power):
        delay = (time-self.end)+(self.end-self.start)/2
        ratio = (self.end-self.start)/(2*delay)
        constant = np.log(delay/scale)
        if power == 1:
            coefficients = -ratio**_DEGREE[1:]/_DEGREE[1:]
        else:
            coefficients = 2*(_HARMONIC[1:]-constant)*ratio**_DEGREE[1:]/_DEGREE[1:]
        return np.tensordot(np.r_[constant**power,coefficients],self.moments,axes=1)

    def distant_inverse(self,time,source_index,power):
        half = (self.end-self.start)/2
        delay = (time-self.end)+half
        ratio = half/delay
        if ratio>1/8:
            return None
        coefficients = _INVERSE_WEIGHTS[power]*ratio**_DEGREE/delay**power
        return np.tensordot(coefficients,self.moments[:,source_index],axes=1)


@dataclass(frozen=True)
class CompressedPanelHistory:
    """Keep every accepted panel and accelerate only the smooth far kernel."""

    start: float
    source_shape: tuple
    identity: str
    incoming: object
    roots: tuple = ()

    @classmethod
    def from_history(cls, history):
        return cls(history.start,history.source_shape,history.identity,history.incoming).append(*history.panels)

    def append(self,*panels):
        end = self.end
        updated = list(self.roots)
        for panel in panels:
            if panel.start != end or panel.coefficients.shape[1:] != self.source_shape or len(panel.coefficients)>3:
                raise ValueError("Compression requires adjacent quadratic-or-lower source panels.")
            end = panel.end
            block = _Block.leaf(panel)
            while updated and updated[-1].count == block.count:
                block = _Block.merge(updated.pop(),block)
            updated.append(block)
        return CompressedPanelHistory(self.start,self.source_shape,self.identity,self.incoming,tuple(updated))

    @property
    def end(self):
        return self.roots[-1].end if self.roots else self.start

    @cached_property
    def panels(self):
        stack = list(reversed(self.roots))
        leaves = []
        while stack:
            node = stack.pop()
            if node.count == 1:
                leaves.append(node.members[0])
            else:
                stack.extend(reversed(node.members))
        return tuple(leaves)

    @cached_property
    def _parts_cache(self):
        return {}

    def integral(self,time,scale,power=1,trial=None):
        if (not np.isfinite([time,scale]).all() or scale<=0 or power not in (0,1,2)
                or not self.start <= time <= (self.end if trial is None else trial.end)):
            raise ValueError("Supply a covered cut, positive scale and log power zero to two.")
        if trial is not None and (trial.start!=self.end or trial.coefficients.shape[1:]!=self.source_shape):
            raise ValueError("The provisional source panel must join the committed history.")
        result = self._committed_integral(time, scale, power)
        if trial is not None and trial.start<time:
            result += trial.integral(time,scale,power)
        return result

    def _committed_integral(self, time, scale, power):
        """Integrate the committed prefix up to min(time, end), including for a trial cut."""
        result = np.zeros(self.source_shape)
        if self.incoming is not None:
            result = result+self.incoming.integral(time,scale,power)
        stack = list(reversed(self.roots))
        while stack:
            node = stack.pop()
            if node.start >= time:
                continue
            if power == 0 and node.end <= time:
                result += node.moments[0]
            elif node.count == 1:
                result += node.members[0].integral(time,scale,power)
            elif time>=node.end and (node.end-node.start)<=2*_RATIO*((time-node.end)+(node.end-node.start)/2):
                # Taylor remainder of log is at most r**25/(25*(1-r)).
                # The exact log-squared coefficients use 2 H_(n-1)/n.
                result += node.contracted(time,scale,power)
            else:
                stack.extend(reversed(node.members))
        return result

    def value(self,time):
        panel = self.covering_panel(time)
        fraction = (time-panel.start)/(panel.end-panel.start)
        return np.polynomial.polynomial.polyval(fraction,panel.coefficients)

    def covering_panel(self,time):
        if not np.isfinite(time) or not self.start<=time<=self.end or not self.roots:
            raise ValueError("A committed panel must cover the queried time.")
        block = next(node for node in self.roots if time<=node.end)
        while block.count>1:
            left,right = block.members
            block = left if time<=left.end else right
        return block.members[0]

    def integrate_jet(self,time,scale,power,orders):
        return self.integral(time,scale,power)[orders]

    def start_jump(self,order):
        return self.value(self.start)[order]

    def intersecting(self,start,stop):
        """Visit only the accepted leaves that intersect the clipped range."""
        stack = list(reversed(self.roots))
        while stack:
            block = stack.pop()
            if block.end<=start or block.start>=stop:
                continue
            if block.count==1:
                yield block.members[0]
            else:
                stack.extend(reversed(block.members))

    def inverse_four(self,time,stop,source_index,leaf_integral):
        """Integrate the positive-delay inverse quartic up to the far cut."""
        return self.inverse_power(time,stop,source_index,leaf_integral,4)

    def inverse_power(self,time,stop,source_index,leaf_integral,power):
        """Contract a far inverse power or its first derivative kernel."""

        total = np.zeros(self.source_shape[1:])
        stack = list(reversed(self.roots))
        while stack:
            block = stack.pop()
            if block.start>=stop:
                continue
            if block.count==1:
                total += leaf_integral(block.members[0],time,stop,source_index)
                continue
            part = block.distant_inverse(time,source_index,power) if block.end<=stop else None
            if part is not None:
                total += part
            else:
                stack.extend(reversed(block.members))
        return total
