from asdfspline import PiecewiseMonotoneCubicSpline
import numpy as np
import pytest


# TODO: separate test functions?
def test_values_errors():
    with pytest.raises(ValueError, match='at least two values'):
        PiecewiseMonotoneCubicSpline(0)
    with pytest.raises(ValueError, match='at least two values'):
        PiecewiseMonotoneCubicSpline([0])
    with pytest.raises(ValueError, match='grid.*same as.*values'):
        PiecewiseMonotoneCubicSpline([0, 1], grid=[0], closed=False)
    with pytest.raises(ValueError, match='grid.*one more than.*values'):
        PiecewiseMonotoneCubicSpline([0, 1], grid=[0, 1], closed=True)


def test_grid_errors():
    with pytest.raises(ValueError, match='NaN.*not allowed in grid'):
        PiecewiseMonotoneCubicSpline([0, 1], grid=[0, np.NaN])
    with pytest.raises(ValueError, match='grid.*must be strictly ascending'):
        PiecewiseMonotoneCubicSpline([0, 1], grid=[0, 0])


def test_slopes_errors():
    with pytest.raises(ValueError, match='slopes.*same as.*values'):
        PiecewiseMonotoneCubicSpline([0, 1], slopes=[0])
    with pytest.raises(ValueError, match='slope.*index 0.*too steep'):
        PiecewiseMonotoneCubicSpline([0, 1], slopes=[4, 0])
    with pytest.raises(ValueError, match='slope.*index 0.*wrong sign'):
        PiecewiseMonotoneCubicSpline([0, 1], slopes=[-1, 0])


def test_python_errors():
    with pytest.raises(ValueError, match='one-dimensional'):
        PiecewiseMonotoneCubicSpline([[0]])
    with pytest.raises(ValueError, match='grid.*one-dimensional'):
        PiecewiseMonotoneCubicSpline([0, 1], grid=[[0]])
    with pytest.raises(ValueError, match='slopes.*one-dimensional'):
        PiecewiseMonotoneCubicSpline([0, 1], slopes=[[0]])
