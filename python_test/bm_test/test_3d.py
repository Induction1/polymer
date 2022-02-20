import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from bm.brownian_motion import calculate_f, from_bond_vectors_of_a_chain_to_padded_bead_vectors, \
    generate_3d_random_bm_vectors, generate_random_3d_unit_vectors_of_n_chains, \
    from_bond_vectors_of_n_chains_to_padded_bead_vectors, vector_magnitude


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
        result = [self.single_dna_chain_force_extension(steps=steps) for steps in [10_000, 1000_000, 5000_000]]

        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        # Using Python unpack syntax. Could also use cum_sum[0], cum_sum[1], cum_sum[2]
        for beads_0_to_n_plus_1 in result:
            ax.plot(*beads_0_to_n_plus_1, label='DNA chain')
            ax.legend()
            ax.scatter(*beads_0_to_n_plus_1)
        plt.show()

    def single_dna_chain_force_extension(self, steps=400_000, delta_t=1.0 / 500_000, n=50, k=500, p_magnitude=100):
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
        :param steps:
        :param delta_t:
        :param n:
        :param k:
        :param p_magnitude:
        :return:
        """

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
            n_plus_2_bonds = np.diff(n_plus_3_beads, axis=1)

            n_plus_1_bm = generate_3d_random_bm_vectors(n + 1)
            beads_0_to_n_plus_1 = (n_plus_3_beads[:, 1:n + 2]
                                   - spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, :n + 1])
                                   + spring_switch * delta_t * k * calculate_f(n_plus_2_bonds[:, 1:n + 2])
                                   + p_switch * delta_t * p
                                   + bm_switch * np.sqrt(2 * delta_t) * n_plus_1_bm
                                   )
            bond_vectors = np.diff(beads_0_to_n_plus_1, axis=1)

        return beads_0_to_n_plus_1

    def test_force_extension_of_dna_multi_chains(self):
        num_of_chains = 50
        steps = 500000
        delta_t = 0.0001

        k = 1000
        n_array = [10, 20, 50, 100, 150, 200]

        p_array = [0, 2, 5, 10, 20, 50, 100, 200, 500]
        print('len,k,p,num_of_chains,steps,delta_t,time,chain_e2e_avg_distance')

        for p_magnitude in p_array:
            for n in n_array:
                spring_switch = 1  # Turn on/off (1/0) the effect of the spring
                p_switch = 1  # Turn on/off (1/0) the effect of the force p
                bm_switch = 1  # Turn on/off (1/0) the effect of Brownian Motion
                bm_factor = 1  # Be able to adjust the magnitude of the random kicks

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

                bead_0 = ensemble_chain[:, 0]
                bead_n = ensemble_chain[:, -1]


                print(f'{n},{k},{p_magnitude},{num_of_chains},{steps},{delta_t},{steps * delta_t},{vector_magnitude(bead_n - bead_0)}')

        # mpl.rcParams['legend.fontsize'] = 10
        # fig = plt.figure()
        # ax = fig.gca(projection='3d')
        # # Using Python unpack syntax. Could also use array[0], array[1], array[2]
        # ax.plot(*ensemble_chain, label='DNA chain')
        # ax.plot(*beads_0_to_n_plus_1[0], label='DNA chain')
        #
        # ax.legend()
        # ax.scatter(*ensemble_chain)
        # ax.scatter(*beads_0_to_n_plus_1[0])
        # plt.show()

    def test_vector_magnitude(self):
        v = np.array([3, 4])
        print(vector_magnitude(v))
