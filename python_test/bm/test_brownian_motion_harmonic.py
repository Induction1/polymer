import unittest

import numpy as np
import matplotlib.pyplot as plt


class TestBMHarmonic(unittest.TestCase):
    """
    Brownian Motion in a harmonic potential
    """
    def test_harmonic(self):
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
        """

        T = 12
        steps = 100 * T
        times = np.linspace(0, T, steps)

        msd = 3 * (1 - np.exp(-2 * times))

        plt.figure(figsize=(9, 6))
        plt.plot(times, msd)
        plt.show()

    def test(self):
        dimensions = 1000
        T = 12
        steps = 100 * T
        times = np.linspace(0, T, steps)

        dt = times[1] - times[0]

        dB = np.sqrt(2 * dt) * np.random.normal(size=(steps, dimensions))
        B = np.zeros(shape=(steps, dimensions))

        r_curr = np.zeros(dimensions)
        for i in range(1, steps):
            B[i] = r_curr * (1.0 - dt) + dB[i]
            r_curr = B[i]

        B = B * B
        msd_3d = 3 * np.sum(B, axis=1) / dimensions
        print([b.shape for b in [dB, B, msd_3d]])
        plt.figure(figsize=(9, 6))
        plt.plot(times, msd_3d)
        plt.show()
