from _asdfspline import ffi as _ffi, lib as _lib
import numpy as _np


class _FromPtr:

    def __init__(self, ptr):
        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
        self._ptr = ptr


class AsdfSpline(_FromPtr):

    def __init__(self, data):
        positions = []
        times = []
        speeds = []
        tcb = []
        closed = False
        for vertex in data:
            vertex = dict(vertex)  # Make a copy
            position = vertex.pop('position', None)
            if position is None:
                raise ValueError('every vertex must have a position')
            time = vertex.pop('time', _np.nan)
            if not _np.isscalar(time):
                raise TypeError('time values must be scalars')
            times.append(time)
            if position == 'closed':
                if vertex:
                    raise ValueError(
                        'when position is "closed", only time is allowed')
                closed = True
                break
            position = _np.asarray(position, dtype='float32')
            if position.ndim != 1:
                raise ValueError('positions must be one-dimensional')
            if len(position) == 2:
                position = _np.append(position, 0)
            if len(position) != 3:
                raise ValueError('positions must be 2- or 3-dimensional')
            positions.append(position)
            speed = vertex.pop('speed', _np.nan)
            if not _np.isscalar(speed):
                raise TypeError('speed values must be scalars')
            speeds.append(speed)
            tension = vertex.pop('tension', 0)
            if not _np.isscalar(tension):
                raise TypeError('tension values must be scalars')
            continuity = vertex.pop('continuity', 0)
            if not _np.isscalar(continuity):
                raise TypeError('continuity values must be scalars')
            bias = vertex.pop('bias', 0)
            if not _np.isscalar(bias):
                raise TypeError('bias values must be scalars')
            tcb.append((tension, continuity, bias))
            if vertex:
                raise ValueError('invalid key(s): {}'.format(set(vertex)))
        if len(positions) < 2:
            raise ValueError('at least two vertices are needed')
        if not closed:
            assert len(tcb) >= 2
            if tcb.pop(0) != (0, 0, 0):
                raise ValueError(
                    'first vertex cannot have tension/continuity/bias '
                    '(except for closed curves')
            if tcb.pop(-1) != (0, 0, 0):
                raise ValueError(
                    'last vertex cannot have tension/continuity/bias '
                    '(except for closed curves)')
        positions, positions_ptr = _make_buffer(3, positions)
        times, times_ptr = _make_buffer(1, times)
        speeds, speeds_ptr = _make_buffer(1, speeds)
        tcb = _np.ascontiguousarray(tcb, dtype='float32')
        tcb_ptr = _ffi.from_buffer('float[]', tcb)
        ptr = _ffi.gc(
            _lib.asdf_asdfposspline(
                positions_ptr, len(positions),
                times_ptr, len(times),
                speeds_ptr, len(speeds),
                tcb_ptr, len(tcb),
                closed,
                ),
            _lib.asdf_asdfposspline_free)
        super().__init__(ptr)

    def evaluate(self, t):
        return _evaluate(t, (3,), _lib.asdf_asdfposspline_evaluate, self._ptr)

    @property
    def grid(self):
        return _grid(_lib.asdf_asdfposspline_grid, self._ptr)


def _evaluate(t, extra_dim, func, ptr):
    t = _np.ascontiguousarray(t, dtype='float32')
    t_ptr = _ffi.from_buffer('float[]', t)
    output = _np.empty(t.shape + extra_dim, dtype='float32')
    output_ptr = _ffi.from_buffer('float[]', output)
    func(ptr, t_ptr, t.size, output_ptr)
    return output


def _grid(func, ptr):
    grid_ptr = _ffi.new('float**')
    grid_len = func(ptr, grid_ptr)
    if grid_len == 0:
        raise RuntimeError(_ffi.string(_lib.asdf_last_error()).decode())
    buffer = _ffi.buffer(
        grid_ptr[0], grid_len * _np.dtype('float32').itemsize)
    array = _np.frombuffer(buffer, dtype='float32')
    # NB: Writing to this array would be Undefined Behavior
    array.flags.writeable = False
    return array


def _make_buffer(dim, numbers, name=None):
    # NB: If name is None, this is supposed to never fail
    numbers = _np.ascontiguousarray(numbers, dtype='float32')
    if dim > 1:
        if numbers.ndim != 2:
            raise ValueError(
                name + ' must be two-dimensional (list of coordinate pairs)')
        if numbers.shape[1] != dim:
            raise ValueError(
                '{} must be a list of {}D coordinate pairs'.format(name, dim))
    elif dim == 1:
        if numbers.ndim != 1:
            raise ValueError(name + ' must be one-dimensional')
    else:
        assert False
    numbers_ptr = _ffi.from_buffer('float[]', numbers)
    return numbers, numbers_ptr


def _check_tcb(tcb, closed, N):
    if tcb is None:
        tcb = 0, 0, 0
    tcb = _np.ascontiguousarray(tcb, dtype='float32')
    if tcb.ndim == 1:
        if not closed:
            N = max(N - 2, 0)
        tcb = _np.tile(tcb, (N, 1))
    if tcb.ndim != 2:
        raise ValueError(
            'TCB values must be two-dimensional (list of triples)')
    if tcb.shape[1] != 3:
        raise ValueError('TCB values must be a list of triples')
    tcb_ptr = _ffi.from_buffer('float[]', tcb)
    return tcb, tcb_ptr


class _CubicCurve3(_FromPtr):

    def evaluate(self, t):
        return _evaluate(t, (3,), _lib.asdf_cubiccurve3_evaluate, self._ptr)

    @property
    def grid(self):
        return _grid(_lib.asdf_cubiccurve3_grid, self._ptr)


class CentripetalKochanekBartelsSpline3(_CubicCurve3):

    def __init__(self, positions, *, tcb=None, closed=False):
        positions, positions_ptr = _make_buffer(3, positions, 'positions')
        tcb, tcb_ptr = _check_tcb(tcb, closed, len(positions))
        ptr = _ffi.gc(
            _lib.asdf_centripetalkochanekbartelsspline3(
                positions_ptr, len(positions),
                tcb_ptr, len(tcb),
                closed,
            ),
            _lib.asdf_cubiccurve3_free)
        super().__init__(ptr)


class _CubicCurve2(_FromPtr):

    def evaluate(self, t):
        return _evaluate(t, (2,), _lib.asdf_cubiccurve2_evaluate, self._ptr)

    @property
    def grid(self):
        return _grid(_lib.asdf_cubiccurve2_grid, self._ptr)


class CentripetalKochanekBartelsSpline2(_CubicCurve2):

    def __init__(self, positions, *, tcb=None, closed=False):
        positions, positions_ptr = _make_buffer(2, positions, 'positions')
        tcb, tcb_ptr = _check_tcb(tcb, closed, len(positions))
        ptr = _ffi.gc(
            _lib.asdf_centripetalkochanekbartelsspline2(
                positions_ptr, len(positions),
                tcb_ptr, len(tcb),
                closed,
            ),
            _lib.asdf_cubiccurve2_free)
        super().__init__(ptr)


class _CubicCurve1(_FromPtr):

    def evaluate(self, t):
        return _evaluate(t, (), _lib.asdf_cubiccurve1_evaluate, self._ptr)

    @property
    def grid(self):
        return _grid(_lib.asdf_cubiccurve1_grid, self._ptr)


class ShapePreservingCubicSpline(_CubicCurve1):

    def __init__(self, values, *, slopes=None, grid=None, closed=False):
        values, values_ptr = _make_buffer(1, values, 'values')
        if grid is None:
            grid = _np.arange(len(values) + closed, dtype='float32')
        grid, grid_ptr = _make_buffer(1, grid, 'grid')
        if slopes is None:
            ptr = _ffi.gc(
                _lib.asdf_shapepreservingcubicspline(
                    values_ptr, len(values),
                    grid_ptr, len(grid),
                    closed,
                ),
                _lib.asdf_cubiccurve1_free)
        else:
            slopes, slopes_ptr = _make_buffer(1, slopes, 'slopes')
            ptr = _ffi.gc(
                _lib.asdf_shapepreservingcubicspline_with_slopes(
                    values_ptr, len(values),
                    slopes_ptr, len(slopes),
                    grid_ptr, len(grid),
                    closed,
                ),
                _lib.asdf_cubiccurve1_free)
        super().__init__(ptr)


class MonotoneCubicSpline(_CubicCurve1):

    def __init__(self, values, *, slopes=None, grid=None):
        values, values_ptr = _make_buffer(1, values, 'values')
        if grid is None:
            grid = _np.arange(len(values), dtype='float32')
        grid, grid_ptr = _make_buffer(1, grid, 'grid')
        if slopes is None:
            ptr = _ffi.gc(
                _lib.asdf_monotonecubic(
                    values_ptr, len(values),
                    grid_ptr, len(grid),
                ),
                _lib.asdf_monotonecubic_free)
        else:
            slopes, slopes_ptr = _make_buffer(1, slopes, 'slopes')
            ptr = _ffi.gc(
                _lib.asdf_monotonecubic_with_slopes(
                    values_ptr, len(values),
                    slopes_ptr, len(slopes),
                    grid_ptr, len(grid),
                ),
                _lib.asdf_monotonecubic_free)
        super().__init__(ptr)  # Just for error-handling
        self._monotone_ptr = ptr
        super().__init__(_lib.asdf_monotonecubic_inner(ptr))

    def get_time(self, values):
        values = _np.ascontiguousarray(values, dtype='float32')
        output = _np.empty_like(values)
        _lib.asdf_monotonecubic_get_time(
            self._monotone_ptr,
            _ffi.from_buffer('float[]', values),
            values.size,
            _ffi.from_buffer('float[]', output))
        return output
