"""Check the exact STF storage used by the high-rank hereditary histories."""

import unittest

import numpy as np

from eccpn.compact_moments import CompactMomentHistory, pack_stf, unpack_stf
from eccpn.history import PanelHistory, SourcePanel
from eccpn.projection import stf_of_vectors


class CompactMoments(unittest.TestCase):
    def test_unpacked_values_and_logs_retain_all_stf_components(self):
        random = np.random.default_rng(6147)
        times = (-1.4,0.2,1.3,1.9)
        for rank in range(3,8):
            with self.subTest(rank=rank):
                tensors = stf_of_vectors([random.uniform(-1,1,(3,2,3)) for _ in range(rank)])
                self.assertEqual(pack_stf(tensors,rank).shape,(3,2,2*rank+1))
                np.testing.assert_allclose(unpack_stf(pack_stf(tensors,rank),rank),tensors,atol=1e-14,rtol=1e-13)

                panels = [
                    SourcePanel(a,b,(1+k)*tensors)
                    for k,(a,b) in enumerate(zip(times[:-1],times[1:]))
                ]
                direct = PanelHistory(times[0],tensors.shape[1:],"full-stf",tuple(panels))
                packed = tuple(SourcePanel(item.start,item.end,pack_stf(item.coefficients,rank))
                               for item in panels)
                base = PanelHistory(times[0],packed[0].coefficients.shape[1:],"packed-stf",packed[:2])
                suffix = PanelHistory(times[-2],base.source_shape,base.identity,packed[2:])
                compact = CompactMomentHistory(base,rank).extend(CompactMomentHistory(suffix,rank))
                np.testing.assert_allclose(compact.value(0.6),direct.value(0.6),rtol=1e-13,atol=1e-14)
                for order in (0,1,2):
                    for cut in (0.6,1.6,1.9):
                        packed_logs = compact.integrate_jet(cut,0.48,order,slice(None))
                        direct_logs = direct.integrate_jet(cut,0.48,order,slice(None))
                        np.testing.assert_allclose(packed_logs,direct_logs,rtol=1e-13,atol=1e-14)


if __name__ == "__main__":
    unittest.main()
