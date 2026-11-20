"""Canonical and orbital histories must depend only on the physical inputs."""
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from eccpn import Accuracy, Binary, InitialState, evolve
from eccpn.hereditary import AcceptedHistory
from eccpn.history import FiniteOrbitalPast
from eccpn.incoming import IncomingOrbit, ParabolicOrbit


class ParallelHistory(unittest.TestCase):
    binary = Binary(.68, .32, c=1.4)
    accuracy = Accuracy(max_step=.13, rtol=2e-8, incoming_waveform_nodes=8)

    def finite(self):
        initial = InitialState((20., 0., 0.), (.073, .251, .009))
        precursor = evolve(self.binary, initial, None, (-.4, 0.), Accuracy("newtonian"))
        initial = InitialState(*precursor.state(0.))
        angular_momentum = self.binary.mu * np.cross(initial.position, initial.velocity)
        past = FiniteOrbitalPast(
            (-.4, -.17, 0.), precursor.state, {0: 1., 2: -.002, 4: -.00003},
            {0: angular_momentum, 2: .02 * angular_momentum}, 1., "parallel-finite-prefix",
        )
        zero = (0., 0., 0.)
        return initial, past.initial_data(self.binary, zero, zero, zero)

    def compare(self, observed, reference):
        np.testing.assert_array_equal(observed._solution.boundaries, reference._solution.boundaries)
        np.testing.assert_array_equal(observed._solution.dense, reference._solution.dense)
        self.assertEqual(observed._solution.rejected_steps, reference._solution.rejected_steps)
        actual, expected = observed.checkpoint, reference.checkpoint
        np.testing.assert_array_equal(actual.state, expected.state)
        self.assertEqual(actual.next_step, expected.next_step)
        self.assertEqual(actual.configuration, expected.configuration)
        self.assertEqual(len(actual.history.panels), len(expected.history.panels))
        for left, right in zip(actual.history.panels, expected.history.panels):
            self.assertEqual((left.start, left.end), (right.start, right.end))
            np.testing.assert_array_equal(left.coefficients, right.coefficients)
        self.assertEqual(actual.canonical_history.to_record(), expected.canonical_history.to_record())

    def test_physical_incoming_histories(self):
        orbits = (IncomingOrbit((.02, .17, .03), (23., 0., 0.), -1300.),
                  ParabolicOrbit(24., (1., 0., 0.), (0., 1., 0.), -1300.))
        with ProcessPoolExecutor(2, mp_context=multiprocessing.get_context("spawn")) as pool:
            for orbit in orbits:
                with self.subTest(orbit=type(orbit).__name__):
                    initial, past = orbit.prepare(self.binary)
                    span = (orbit.cutoff, orbit.cutoff + .56)
                    reference = evolve(self.binary, initial, past, span, self.accuracy)
                    actual = evolve(self.binary, initial, past, span, self.accuracy, executor=pool)
                    self.compare(actual, reference)

    def test_finite_past_and_checkpoint_across_workers(self):
        initial, past = self.finite()
        first_span, second_span = (0., .47), (.47, .94)
        reference = evolve(self.binary, initial, past, first_span, self.accuracy)
        with ProcessPoolExecutor(2, mp_context=multiprocessing.get_context("spawn")) as pool:
            first = evolve(self.binary, initial, past, first_span, self.accuracy, executor=pool, parallel_orbit=True)
        self.compare(first, reference)
        with TemporaryDirectory() as directory:
            path = Path(directory) / "orbit.json"
            first.checkpoint.save(path)
            restored = AcceptedHistory.load(path, identity=past.quadrupole.identity)
        state = InitialState(*reference.state(first_span[1]))
        expected = evolve(self.binary, state, reference.checkpoint, second_span, self.accuracy)
        with ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn")) as pool:
            actual = evolve(self.binary, state, restored, second_span, self.accuracy, executor=pool, parallel_orbit=True)
        self.compare(actual, expected)


if __name__ == "__main__":
    unittest.main()
