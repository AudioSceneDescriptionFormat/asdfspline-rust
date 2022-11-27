import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


def plot_1d(s, ax=None):
    """Plot a one-dimensional spline."""
    if ax is None:
        ax = plt.gca()
    times = np.linspace(s.grid[0], s.grid[-1], 200, endpoint=True)
    ax.plot(times, s.evaluate(times))
    ax.scatter(s.grid, s.evaluate(s.grid), marker='x', c='black')


def plot_2d(s, dots_per_second=10, ax=None):
    """Plot a two-dimensional spline (or a 3D-spline with all zero y-values)."""
    total_duration = s.grid[-1] - s.grid[0]
    times = s.grid[0] + np.arange(int(total_duration * dots_per_second) + 1) / dots_per_second
    if ax is None:
        ax = plt.gca()
    data = s.evaluate(times)
    if data.shape[1] == 3:
        if np.any(data[:, 2] != 0):
            raise ValueError('z values must be zero')
        data = data[:, :2]
    ax.plot(*data.T, '.')
    ax.scatter(*s.evaluate(s.grid)[:, :2].T, marker='x', c='black')
    ax.axis('equal')


def plot_3d(s, dots_per_second=10, ax=None):
    total_duration = s.grid[-1] - s.grid[0]
    times = s.grid[0] + np.arange(int(total_duration * dots_per_second) + 1) / dots_per_second
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    ax.plot(*s.evaluate(times).T, '.')
    # https://github.com/matplotlib/matplotlib/issues/18020
    #ax.scatter(*s.evaluate(s.grid).T, marker='x', c='black')
    ax.scatter(*s.evaluate(s.grid).T, c='black')
    set_3d_axes_equal(ax)


def set_3d_axes_equal(ax):
    # https://stackoverflow.com/a/50664367/
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])
