import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from bm.brownian_motion import calculate_f, from_bond_vectors_of_a_chain_to_padded_bead_vectors, \
    generate_3d_random_bm_vectors


class ThreeDimTestCase(unittest.TestCase):

    def test_3d_random_walk(self):
        """
        Demonstrate 3D random walk
        """
        n = 50

        # Generate spherically random vectors with step length 1
        data = generate_3d_random_bm_vectors(n)
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

        all_n_plus_3_beads = from_bond_vectors_of_a_chain_to_padded_bead_vectors(bond_vectors)
        self.assertTrue(
            np.array_equal(
                np.array([[-1., 0., 1., 2.],
                          [0., 0., 1., 1.],
                          [0., 0., 1., 1.]]),
                all_n_plus_3_beads)
        )

    def test_force_extension_of_dna(self):
        steps = 25000
        delta_t = 0.05
        n = 30
        k = 0.5
        p_magnitude = 2

        spring_switch = 1  # Turn on/off (1/0) the effect of the spring
        p_switch = 1  # Turn on/off (1/0) the effect of the force p
        bm_switch = 1  # Turn on/off (1/0) the effect of Brownian Motion

        bond_vectors = generate_3d_random_bm_vectors(n)
        p = np.zeros(3 * (n + 1)).reshape(3, -1)
        p[0, 0] = -p_magnitude
        p[0, n] = p_magnitude

        for i in range(steps):
            n_plus_3_beads = from_bond_vectors_of_a_chain_to_padded_bead_vectors(bond_vectors)
            n_plus_2_bonds = np.diff(n_plus_3_beads, axis=1)

            n_plus_1_bm = generate_3d_random_bm_vectors(n + 1)
            beads_0_to_n_plus_1 = (n_plus_3_beads[:, 1:n + 2]
                                   - spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, :n + 1])
                                   + spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, 1:n + 2])
                                   + p_switch * delta_t * p
                                   + bm_switch * np.sqrt(2 * delta_t) * n_plus_1_bm
                                   )
            bond_vectors = np.diff(beads_0_to_n_plus_1)

        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Using Python unpack syntax. Could also use cum_sum[0], cum_sum[1], cum_sum[2]
        ax.plot(*beads_0_to_n_plus_1, label='DNA chain')
        ax.legend()
        ax.scatter(*beads_0_to_n_plus_1)
        plt.show()
