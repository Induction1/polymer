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


def generate_random_3d_unit_vectors(num_of_vectors):
    """
    :param num_of_vectors: The number of unit random 3D vectors to produce.
    :return: num_of_vectors randomly and uniformly distributed 3D unit vectors.
    """
    angles = np.random.uniform(0.0, 1.0, 2 * num_of_vectors).reshape((2, -1))
    x = np.cos(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    y = np.sin(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    z = np.cos(angles[0] * 2 * np.pi)
    vectors_3d = np.concatenate((x, y, z), axis=0)
    return vectors_3d.reshape(3, -1)


def generate_random_3d_unit_vectors_of_n_chains(num_of_chains, num_of_vectors):
    """

    :param num_of_chains:
    :param num_of_vectors: The number of unit random 3D vectors to produce.
    :return:
    """

    angles = np.random.uniform(0.0, 1.0, 2 * num_of_chains * num_of_vectors).reshape((2, -1))
    x = np.cos(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    y = np.sin(angles[1] * 2 * np.pi) * np.sin(angles[0] * 2 * np.pi)
    z = np.cos(angles[0] * 2 * np.pi)
    vectors_3d = np.concatenate((x, y, z), axis=0)
    return vectors_3d.reshape(num_of_chains, 3, -1)


def vector_magnitude(v):
    return np.sqrt(np.sum(v * v, axis=0))


def calculate_f(v):
    """
    See definition in problem II.1
    f(vector) = vector * (1 - 1/|vector|)
    :param v:
    :return:
    """
    try:
        magnitude = np.sqrt(np.sum(v * v, axis=-2))
        magnitude = np.expand_dims(magnitude, axis=-2)
        return v - v / magnitude
    except:
        return np.zeros(v.shape)


def from_bond_vectors_of_a_chain_to_padded_bead_vectors(bond_vectors_of_a_chain):
    """
    Given n bond_vectors which links n+1 bead, assuming bead_0 is at original, produce the
    coordinates all the n+1 beads plus the padded the before and after bead.
    The before bead is 1 unit to the left of the origin and the after bead is one unit to the
    right of bean_n_plus_1.
    As a result, the coordinates n+3 beads are returned.
    :param bond_vectors_of_a_chain:
    :return:
    """

    dim = bond_vectors_of_a_chain.shape[0]

    bead_neg_1 = np.zeros((dim, 1))
    bead_neg_1[0, 0] = -1

    bead_0 = np.zeros((dim, 1))

    bead_1_to_n_plus_1 = np.cumsum(bond_vectors_of_a_chain, axis=1)

    bead_n_plus_2 = np.array(bead_1_to_n_plus_1[:, -1]).reshape(dim, 1)
    bead_n_plus_2[0, 0] = bead_n_plus_2[0, 0] + 1

    all_beads = np.concatenate((bead_neg_1, bead_0, bead_1_to_n_plus_1, bead_n_plus_2), axis=1)
    return all_beads


def from_bond_vectors_of_n_chains_to_padded_bead_vectors(bond_vectors_of_n_chains):
    """

    :param bond_vectors_of_n_chains: Represents the bond vectors of n chains. Each has vectors of d dimension.
    Each vector has l length. The shape is therefore [n, d, l]
    :return: The padded chains. Each chain is padded with a before-bead that is 1 unit to the left of bead_0 and an
    after-bead that is 1 unit to the right of the last bead bead_-1. Overall, each chain has n+3 beads (Note, the input
    has n vectors, which connect n+1 beads. After padding with before/after beads, the total number of beads of each
    chain becomes n+3)
    """
    num_of_chains = bond_vectors_of_n_chains.shape[0]
    dim = bond_vectors_of_n_chains.shape[1]

    bead_neg_1 = np.zeros((num_of_chains, dim, 1))
    bead_neg_1[:, 0, 0] = -1

    bead_0 = np.zeros((num_of_chains, dim, 1))

    bead_1_to_n_plus_1 = np.cumsum(bond_vectors_of_n_chains, axis=2)

    bead_n_plus_2 = np.array(bead_1_to_n_plus_1[:, :, -1]).reshape(num_of_chains, dim, 1)
    bead_n_plus_2[:, 0, 0] = bead_n_plus_2[:, 0, 0] + 1

    all_chains = np.concatenate((bead_neg_1, bead_0, bead_1_to_n_plus_1, bead_n_plus_2), axis=2)
    return all_chains


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


def generate_simple_3d_unit_vectors(n):
    return np.concatenate((np.ones((n, 1)), np.zeros((n, 2))), axis=1)
