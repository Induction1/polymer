from math import pi, sin, cos
import random
import numpy as np


def generate_random_3d_unit_vector():
    phi = 2 * pi * random.uniform(0, 1)
    theta = 2 * pi * random.uniform(0, 1)
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return x, y, z


def generate_random_3d_unit_vectors(n):
    angles = np.random.uniform(0.0, 1.0, 2 * n).reshape((2, -1))
    x = np.cos(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    y = np.sin(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    z = np.cos(angles[0] * 2 * np.pi)
    vectors_3d = np.concatenate((x, y, z), axis=0)
    return vectors_3d.reshape(3, -1)
