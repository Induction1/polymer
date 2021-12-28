import unittest
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.ticker import PercentFormatter, FuncFormatter


class MyTestCase(unittest.TestCase):
    def test_bm(self):
        """
        Tests Brownian Motion
        :return:
        """

        n = 100
        T = 1.
        d = 100
        times = np.linspace(0, T, n)
        dt = times[1] - times[0]

        # Bt2 - Bt1 ~ Normal with mean 0 and variance t2-t1

        dB = np.sqrt(dt) * np.random.normal(size=(n-1, d))
        BO = np.zeros(shape=(1,d))
        B = np.concatenate((BO, np.cumsum(dB, axis=0)), axis=0)
        print(times.shape)
        print(B.shape)
        plt.plot(times, B)
        plt.show()

    def test_msd(self):
        """
        Tests the mean square distance of a Brownian motion
        :return:
        """

        n = 1000
        T = 1.
        d = 1000
        times = np.linspace(0, T, n)
        dt = times[1]-times[0]

        # Bt2 - Bt1 ~ normal with mean 0 and variance t2 - t1
        dB = np.sqrt(dt) * np.random.normal(size=(n-1,d))
        BO = np.zeros(shape=(1,d))
        B = np.concatenate((BO, np.cumsum(dB, axis=0)), axis=0)

        msd = np.sum(B * B / d, axis=1)
        plt.plot(times, msd)
        plt.show()

    def test_histograph(self):
        n = 1000    # Number of steps
        T = 1.      # Total time in seconds
        d = 10000    # Total number of instances
        times = np.linspace(0, T, n)
        dt = times[1] - times[0]

        # Bt2 - Bt1 ~ normal with mean 0 and variance t2 - t1
        dB = np.sqrt(dt) * np.random.normal(size=(n - 1, d))
        BO = np.zeros(shape=(1, d))
        B = np.concatenate((BO, np.cumsum(dB, axis=0)), axis=0)

        cross_section = B[900, :]

        print(cross_section.shape)

        def fmt(y, position):
            # Ignores the passed in position. This has the effect of scaling the default
            # tick locations.

            s = str(100.0*y / d)

            # The percent symbol needs escaping in LaTeX
            if rcParams['text.usetex'] is True:
                return s + r'$\%$'
            else:
                return s + '%'
        n, bins, patches = plt.hist(cross_section, bins=100)
        plt.gca().yaxis.set_major_formatter(FuncFormatter(fmt))
        plt.show()