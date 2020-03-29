from asdfspline import ShapePreservingCubicSpline
import numpy as np
import pytest


# TODO: separate test functions?
def test_values_errors():
    with pytest.raises(ValueError, match='at least two values'):
        ShapePreservingCubicSpline(0)
    with pytest.raises(ValueError, match='at least two values'):
        ShapePreservingCubicSpline([0])
    with pytest.raises(ValueError, match='grid.*same as.*values'):
        ShapePreservingCubicSpline([0, 1], grid=[0], closed=False)
    with pytest.raises(ValueError, match='grid.*one more than.*values'):
        ShapePreservingCubicSpline([0, 1], grid=[0, 1], closed=True)


def test_grid_errors():
    with pytest.raises(ValueError, match='NaN.*not allowed in grid'):
        ShapePreservingCubicSpline([0, 1], grid=[0, np.NaN])
    with pytest.raises(ValueError, match='grid.*must be strictly ascending'):
        ShapePreservingCubicSpline([0, 1], grid=[0, 0])


def test_slopes_errors():
    with pytest.raises(ValueError, match='slopes.*same as.*values'):
        ShapePreservingCubicSpline([0, 1], slopes=[0])
    with pytest.raises(ValueError, match='slope.*index 0.*too steep'):
        ShapePreservingCubicSpline([0, 1], slopes=[4, 0])
    with pytest.raises(ValueError, match='slope.*index 0.*wrong sign'):
        ShapePreservingCubicSpline([0, 1], slopes=[-1, 0])


def test_python_errors():
    with pytest.raises(ValueError, match='one-dimensional'):
        ShapePreservingCubicSpline([[0]])
    with pytest.raises(ValueError, match='grid.*one-dimensional'):
        ShapePreservingCubicSpline([0, 1], grid=[[0]])
    with pytest.raises(ValueError, match='slopes.*one-dimensional'):
        ShapePreservingCubicSpline([0, 1], slopes=[[0]])
