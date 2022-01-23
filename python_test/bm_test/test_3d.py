import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from bm.brownian_motion import generate_random_3d_unit_vectors, calculate_f, from_bond_vectors_to_padded_bead_vectors


class ThreeDimTestCase(unittest.TestCase):

    def test_3d_random_walk(self):
        """
        Demonstrate 3D random walk
        """
        n = 50

        # Generate spherically random vectors with step length 1
        data = generate_random_3d_unit_vectors(n=n)
        # Create normal random factors for each step
        random_factor = np.random.normal(size=(1, n))
        # Apply the random factors
        data = random_factor * data
        # Construct the random walk path
        cum_sum = np.cumsum(data, axis=1)

        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Using Python unpack syntax. Could also use cum_sum[0], cum_sum[1], cum_sum[2]
        ax.plot(*cum_sum, label='parametriccurve')
        ax.legend()
        ax.scatter(*cum_sum)
        plt.show()

    def test_calculate_f(self):
        a = np.array([[3], [4]])
        self.assertTrue(np.array_equal(np.array([[2.4], [3.2]]), calculate_f(a)))

    def test_from_bond_vectors_to_padded_bead_vectors(self):
        # A single bond vector (1,1,1), which means bead_0 at (0, 0, 0) and bead_1 at (1, 1, 1)
        bond_vectors = np.ones(3).reshape(3, -1)

        all_n_plus_3_beads = from_bond_vectors_to_padded_bead_vectors(bond_vectors)
        self.assertTrue(
            np.array_equal(
                np.array([[-1., 0., 1., 2.],
                          [0., 0., 1., 1.],
                          [0., 0., 1., 1.]]),
                all_n_plus_3_beads)
        )
