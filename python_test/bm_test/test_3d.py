import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from bm.brownian_motion import calculate_f, from_bond_vectors_of_a_chain_to_padded_bead_vectors, \
    generate_3d_random_bm_vectors, generate_random_3d_unit_vectors_of_n_chains, \
    from_bond_vectors_of_n_chains_to_padded_bead_vectors


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
        """
        Tests the forced extension of a single DNA chain.
        We do the following:
        1. Create n 3-d bond vectors
        2. These n vectors implies n+1 beads where bead_0 is at the origin. We then
        pad the n+1 beads with a before and after bead so that the bond length is always 1.
        The length=1 bond will exert no Hooke tension to bead_0 and bead_n_plus_1. The reason to have them is so that
        we can treat all beads 0 to n+1 with the same math operation.
        3. Apply the Euler approximation to get the next round of bead positions.
        4. Caclulate the bond vectors based on the beads and we go back to 2.
        :return:
        """
        steps = 500000
        delta_t = 0.0005
        n = 50
        k = 500
        p_magnitude = 200

        spring_switch = 1  # Turn on/off (1/0) the effect of the spring
        p_switch = 1  # Turn on/off (1/0) the effect of the force p
        bm_switch = 1  # Turn on/off (1/0) the effect of Brownian Motion

        # starting from a straight line will take longer to reach equilibrium
        # bond_vectors = generate_simple_3d_unit_vectors(n)
        bond_vectors = generate_3d_random_bm_vectors(n)

        p = np.zeros(3 * (n + 1)).reshape(3, -1)
        p[0, 0] = -p_magnitude
        p[0, n] = p_magnitude

        for i in range(steps):
            n_plus_3_beads = from_bond_vectors_of_a_chain_to_padded_bead_vectors(bond_vectors)
            n_plus_2_bonds = np.diff(n_plus_3_beads, axis=-1)

            n_plus_1_bm = generate_3d_random_bm_vectors(n + 1)
            beads_0_to_n_plus_1 = (n_plus_3_beads[:, 1:n + 2]
                                   - spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, :n + 1])
                                   + spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, 1:n + 2])
                                   + p_switch * delta_t * p
                                   + bm_switch * np.sqrt(2 * delta_t) * n_plus_1_bm
                                   )
            bond_vectors = np.diff(beads_0_to_n_plus_1, axis=-1)

        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Using Python unpack syntax. Could also use cum_sum[0], cum_sum[1], cum_sum[2]
        ax.plot(*beads_0_to_n_plus_1, label='DNA chain')
        ax.legend()
        ax.scatter(*beads_0_to_n_plus_1)
        plt.show()

    def test_force_extension_of_dna_multi_chains(self):
        num_of_chains = 50
        steps = 500000
        delta_t = 0.0005
        n = 50
        k = 500
        p_magnitude = 200

        spring_switch = 1  # Turn on/off (1/0) the effect of the spring
        p_switch = 1  # Turn on/off (1/0) the effect of the force p
        bm_switch = 1  # Turn on/off (1/0) the effect of Brownian Motion
        bm_factor = 1

        # starting from a straight line will take longer to reach equilibrium
        bond_vectors = generate_random_3d_unit_vectors_of_n_chains(num_of_chains, n)

        p = np.zeros(3 * num_of_chains * (n + 1)).reshape(num_of_chains, 3, -1)
        p[:, 0, 0] = -p_magnitude
        p[:, 0, n] = p_magnitude

        for i in range(steps):
            n_plus_3_beads = from_bond_vectors_of_n_chains_to_padded_bead_vectors(bond_vectors)
            # TODO the below line is too wasty. We should be able to pad bond_vectors directly.
            n_plus_2_bonds = np.diff(n_plus_3_beads, axis=-1)

            n_plus_1_bm = generate_random_3d_unit_vectors_of_n_chains(num_of_chains, n + 1)
            beads_0_to_n_plus_1 = (n_plus_3_beads[:, :, 1:n + 2]
                                   - spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, :, :n + 1])
                                   + spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, :, 1:n + 2])
                                   + p_switch * delta_t * p
                                   + bm_switch * bm_factor * np.sqrt(2 * delta_t) * n_plus_1_bm
                                   )
            bond_vectors = np.diff(beads_0_to_n_plus_1, axis=-1)

        ensemble_chain = np.sum(beads_0_to_n_plus_1, axis=0) / num_of_chains

        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Using Python unpack syntax. Could also use array[0], array[1], array[2]
        ax.plot(*ensemble_chain, label='DNA chain')
        ax.plot(*beads_0_to_n_plus_1[0], label='DNA chain')

        ax.legend()
        ax.scatter(*ensemble_chain)
        ax.scatter(*beads_0_to_n_plus_1[0])
        plt.show()
