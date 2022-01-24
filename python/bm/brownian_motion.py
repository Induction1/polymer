import random
from math import pi, sin, cos

import numpy as np


def generate_random_3d_unit_vector():
    """

    :return: A single random 3D unit vector.
    """
    phi = 2 * pi * random.uniform(0, 1)
    theta = 2 * pi * random.uniform(0, 1)
    x = sin(theta) * cos(phi)
    y = sin(theta) * sin(phi)
    z = cos(theta)
    return x, y, z


def generate_random_3d_unit_vectors(n):
    """
    :param n: The number of unit random 3D vectors to produce.
    :return: n randomly and uniformly distributed 3D unit vectors.
    """
    angles = np.random.uniform(0.0, 1.0, 2 * n).reshape((2, -1))
    x = np.cos(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    y = np.sin(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    z = np.cos(angles[0] * 2 * np.pi)
    vectors_3d = np.concatenate((x, y, z), axis=0)
    return vectors_3d.reshape(3, -1)


def calculate_f(v):
    """
    See definition in problem II.1
    f(vector) = vector * (1 - 1/|vector|)
    :param v:
    :return:
    """
    try:
        return v - v / np.sqrt(np.sum(v * v, axis=0))
    except:
        return np.zeros(v.shape)


def from_bond_vectors_to_padded_bead_vectors(bond_vectors):
    """
    Given n bond_vectors which links n+1 bead, assuming bead_0 is at original, produce the
    coordinates all the n+1 beads plus the padded the before and after bead.
    The before bead is 1 unit to the left of the origin and the after bead is one unit to the
    right of bean_n_plus_1.
    As a result, the coordinates n+3 beads are returned.
    :param bond_vectors:
    :return:
    """

    dim = bond_vectors.shape[0]

    bead_neg_1 = np.zeros(dim).reshape((dim, -1))
    bead_neg_1[0, 0] = -1

    bead_0 = np.zeros(dim).reshape(dim, -1)

    bead_1_to_n_plus_1 = np.cumsum(bond_vectors, axis=1)

    bead_n_plus_2 = np.array(bead_1_to_n_plus_1[:, -1]).reshape(dim, -1)
    bead_n_plus_2[0, 0] = bead_n_plus_2[0, 0] + 1

    all_beads = np.concatenate((bead_neg_1, bead_0, bead_1_to_n_plus_1, bead_n_plus_2), axis=1)
    return all_beads


def generate_3d_random_bm_vectors(n):
    """
    Generate 3D vectors whose lengths are normally distributed and directions are randomly pointed
    in spherical coordinates.
    :param n: The number of 3D vectors
    :return:
    """
    random_vs = generate_random_3d_unit_vectors(n)
    random_factor = np.random.normal(0, 1, n)
    # Apply the random factors
    random_bm_vs = random_factor * random_vs
    return random_bm_vs
