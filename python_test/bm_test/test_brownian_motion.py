import unittest
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.ticker import FuncFormatter
from bm.brownian_motion import generate_random_3d_unit_vector, generate_random_3d_unit_vectors


class TestSimpleBM(unittest.TestCase):
    """
    Simple Brownian Motion tests
    """

    def get_bm_data(self, dimensions=200, steps=100):
        """

        :param dimensions: number of realizations of bm_test
        :param steps: number of steps in each realization
        :return:
        """
        T = 1.
        times = np.linspace(0, T, steps)
        # print(times.shape)
        dt = times[1] - times[0]

        # Bt2 - Bt1 ~ Normal with mean 0 and variance t2-t1

        dB = np.sqrt(dt) * np.random.normal(size=(steps - 1, dimensions))
        B0 = np.zeros(shape=(1, dimensions))
        print(dB.shape)
        print(B0.shape)
        B = np.concatenate((B0, np.cumsum(dB, axis=0)), axis=0)
        return times, B

    def test_bm(self):
        """
        Tests Brownian Motion
        :return:
        """
        times, B = self.get_bm_data()

        print(times.shape)
        print(B.shape)

        return times, B
        # print(times.shape)
        # print(B.shape)
        # plt.plot(times, B)
        # plt.show()

    def test_msd(self):
        """
        Tests the mean square distance of a Brownian motion
        :return:
        """
        d = 100
        times, B = self.get_bm_data(dimensions=d)

        msd = np.sum(B * B / d, axis=1)

        return times, msd
        # plt.plot(times, msd)
        # plt.show()

    def test_histogram(self):
        n = 1000
        d = 1000
        times, B = self.get_bm_data(dimensions=d, steps=n)

        cross_section = B[900, :]

        print(cross_section.shape)

        def fmt(y, position):
            # Ignores the passed in position. This has the effect of scaling the default
            # tick locations.

            s = str(100.0 * y / d)

            # The percent symbol needs escaping in LaTeX
            if rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'

        plt.hist(cross_section, bins=100)
        plt.gca().yaxis.set_major_formatter(FuncFormatter(fmt))
        # plt.show()

    def test_combine_plots(self):
        n = 1000
        d = 1000
        times, B = self.get_bm_data(dimensions=d, steps=n)
        cross_section = B[900, :]
        msd = np.sum(B * B / d, axis=1)

        def fmt(y, position):
            # Ignores the passed in position. This has the effect of scaling the default
            # tick locations.

            s = str(100.0 * y / d)

            # The percent symbol needs escaping in LaTeX
            if rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'

        fig = plt.figure(figsize=(9, 6))
        sub1 = fig.add_subplot(2, 2, (1, 2))
        sub1.plot(times, B)
        sub2 = fig.add_subplot(2, 2, 3)
        sub2.plot(times, msd)
        sub3 = fig.add_subplot(2, 2, 4)
        sub3.hist(cross_section, bins=100)
        plt.gca().yaxis.set_major_formatter(FuncFormatter(fmt))
        plt.suptitle("Simple Brownian Motion Simulation")
        plt.show()

    def test_random_3d_vector(self):
        for i in range(0, 10):
            x, y, z = generate_random_3d_unit_vector()
            print(x, y, z)
            print(x * x + y * y + z * z)

    def test_random_3d_vector_array(self):
        vectors_3d = generate_random_3d_unit_vectors(4)
        print(vectors_3d.reshape(3, -1))
