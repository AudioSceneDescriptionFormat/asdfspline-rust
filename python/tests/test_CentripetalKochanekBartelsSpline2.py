from asdfspline import CentripetalKochanekBartelsSpline2
import pytest


def test_positions_errors():
    with pytest.raises(ValueError, match='at least two positions'):
        CentripetalKochanekBartelsSpline2([[0, 0]])


def test_tcb_errors():
    with pytest.raises(ValueError, match='positions.*two more than.*TCB'):
        CentripetalKochanekBartelsSpline2(
            [[0, 0], [1, 1]], tcb=[[0, 0, 0]], closed=False)
    with pytest.raises(ValueError, match='positions.*the same.*TCB'):
        CentripetalKochanekBartelsSpline2(
            [[0, 0], [1, 1]], tcb=[[0, 0, 0]], closed=True)


def test_python_errors():
    with pytest.raises(ValueError, match='two-dimensional'):
        CentripetalKochanekBartelsSpline2(0)
    with pytest.raises(ValueError, match='two-dimensional'):
        CentripetalKochanekBartelsSpline2([0])
    with pytest.raises(ValueError, match='list of 2D coordinate pairs'):
        CentripetalKochanekBartelsSpline2([[0]])
    with pytest.raises(ValueError, match='TCB.*two-dimensional'):
        CentripetalKochanekBartelsSpline2([[0, 0], [1, 1]], tcb=[[[0, 0, 0]]])
    with pytest.raises(ValueError, match='TCB.*list of triples'):
        CentripetalKochanekBartelsSpline2([[0, 0], [1, 1]], tcb=0)
    with pytest.raises(ValueError, match='TCB.*list of triples'):
        CentripetalKochanekBartelsSpline2([[0, 0], [1, 1]], tcb=[0])
