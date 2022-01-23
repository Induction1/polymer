import unittest
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from bm.brownian_motion import generate_random_3d_unit_vectors


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
