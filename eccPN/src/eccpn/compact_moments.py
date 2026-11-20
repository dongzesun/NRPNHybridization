"""Exact Cartesian STF storage for canonical hereditary panels."""

from dataclasses import dataclass
from functools import cache, cached_property
from itertools import product

import numpy as np

from .history import SourcePanel


@cache
def _basis(rank):
    independent = tuple((x,rank-z-x,z) for z in (0,1) for x in range(rank-z+1))
    selectors = tuple(np.array([
        [i for i,count in enumerate(counts) for _ in range(count)] for counts in independent
    ]).T)
    vectors = dict(zip(independent,np.eye(len(independent),dtype=np.int8)))

    def vector(counts):
        if counts not in vectors:
            x,y,z = counts
            vectors[counts] = -vector((x+2,y,z-2))-vector((x,y+2,z-2))
        return vectors[counts]

    expansion = np.stack([vector(tuple(indices.count(i) for i in range(3)))
                          for indices in product(range(3),repeat=rank)])
    return selectors,expansion


def pack_stf(tensor,rank):
    selectors,_ = _basis(rank)
    return tensor[(...,*selectors)]


def unpack_stf(coordinates,rank):
    _,expansion = _basis(rank)
    expanded = coordinates @ expansion.T
    return expanded.reshape(coordinates.shape[:-1]+(3,)*rank)


@dataclass(frozen=True)
class _ExpandedPanels:
    panels: tuple
    rank: int

    def __len__(self):
        return len(self.panels)

    def __getitem__(self,index):
        panel = self.panels[index]
        return SourcePanel(panel.start,panel.end,unpack_stf(panel.coefficients,self.rank))

    def __iter__(self):
        for panel in self.panels:
            yield SourcePanel(panel.start,panel.end,unpack_stf(panel.coefficients,self.rank))


@dataclass(frozen=True)
class CompactMomentHistory:
    """Preserve all STF information while exposing the normal tensor interface."""

    packed: object
    rank: int

    @property
    def start(self):
        return self.packed.start

    @property
    def end(self):
        return self.packed.end

    @property
    def source_shape(self):
        return self.packed.source_shape[:-1]+(3,)*self.rank

    @cached_property
    def panels(self):
        return _ExpandedPanels(self.packed.panels,self.rank)

    def value(self,time):
        return unpack_stf(self.packed.value(time),self.rank)

    def start_jump(self,order):
        return unpack_stf(self.packed.start_jump(order),self.rank)

    def integrate_jet(self,time,scale,power,orders):
        return unpack_stf(self.packed.integrate_jet(time,scale,power,orders),self.rank)

    def extend(self,other):
        if self.rank != other.rank:
            raise ValueError("Stored STF moment ranks differ.")
        return CompactMomentHistory(self.packed.append(*other.packed.panels),self.rank)
