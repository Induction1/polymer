import unittest

import numpy as np
import matplotlib.pyplot as plt


class TestBMHarmonic(unittest.TestCase):
    """
    Brownian Motion in a harmonic potential
    """

    def test_harmonic(self):
        times, simulated_msd_3d, msd_3d = TestBMHarmonic.simulate_bm_harmonic(seconds=12)

        plt.figure(figsize=(9, 6))
        plt.plot(times, msd_3d)
        plt.plot(times, simulated_msd_3d)

        plt.show()

    @staticmethod
    def simulate_bm_harmonic(dimensions=1000, seconds=12):
        """
        The below math is normalized by
        t_tilda = t / tao
        r_tilda = r / lambda

        where
            tao = zeta / k and lambda = sqrt(k_B * T / k)
        where
            zeta is the friction coefficient
            k is the Hooke's constant
            k_B is the Boltzmann constant
            T is the absolute temperature

        Analytic MSD solution:
        E[r(t_tilda)^2] = 3 * ( 1 - e^(-2 * t_tilda))
        ---

        Step simulation using Taylor expansion:
        r_tilda_(i+1) = r_tilda_i * ( 1 - delta_t_tilda) + sqrt(2 * delta_t_tilda) * b_i
        Where r_n is the displacement of the nth step, and b_i is
        a stochastic variable of standard normal distribution.

        :param dimensions: The number of 1-d simulation to perform
        :param seconds: The length of the simulation
        :return: times, msd_3d, simulated_msd_3d
        """

        T = 1
        steps = 10 * 100 * T
        times = np.linspace(0, T, steps)
        times = seconds * times
        delta_t = times[1] - times[0]
        delta_brownian_motion = np.sqrt(2 * delta_t) * np.random.normal(size=(steps, dimensions))

        r = np.zeros(shape=(steps, dimensions))

        r_curr = np.zeros(dimensions)
        for i in range(1, steps):
            r[i] = r_curr * (1.0 - delta_t) + delta_brownian_motion[i]
            r_curr = r[i]

        r_square = r * r
        simulated_msd_3d = 3 * (np.sum(r_square, axis=1) / dimensions)

        # The analytical msd_3d
        msd_3d = 3 * (1 - np.exp(-2 * times))

        return times, msd_3d, simulated_msd_3d
