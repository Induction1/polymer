import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def show_empty_3d_grid():
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    plt.show()


def show_3d_spirl():
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
    z = np.linspace(-2, 2, 100)
    r = z ** 2 + 1
    x = r * np.sin(theta)
    y = r * np.cos(theta)
    ax.plot(x, y, z, label='parametriccurve')
    ax.legend()
    plt.show()
