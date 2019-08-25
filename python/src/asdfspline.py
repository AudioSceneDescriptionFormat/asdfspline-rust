from _asdfspline import ffi as _ffi, lib as _lib
import numpy as _np


class AsdfSpline:

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
                raise ValueError('Every vertex must have a position')
            time = vertex.pop('time', _np.nan)
            if not _np.isscalar(time):
                raise TypeError('Time values must be scalars')
            times.append(time)
            if position == 'closed':
                if vertex:
                    raise ValueError(
                        'When position is "closed", only time is allowed')
                closed = True
                break
            position = _np.asarray(position, dtype='float32')
            if position.ndim != 1:
                raise ValueError("Positions must be one-dimensional")
            if len(position) == 2:
                position = _np.append(position, 0)
            if len(position) != 3:
                raise ValueError("Positions must be 2- or 3-dimensional")
            positions.append(position)
            speed = vertex.pop('speed', _np.nan)
            if not _np.isscalar(speed):
                raise TypeError('Speed values must be scalars')
            speeds.append(speed)
            tension = vertex.pop('tension', 0)
            if not _np.isscalar(tension):
                raise TypeError('Tension values must be scalars')
            continuity = vertex.pop('continuity', 0)
            if not _np.isscalar(continuity):
                raise TypeError('Continuity values must be scalars')
            bias = vertex.pop('bias', 0)
            if not _np.isscalar(bias):
                raise TypeError('Bias values must be scalars')
            tcb.append((tension, continuity, bias))
            if vertex:
                raise ValueError('Invalid key(s): {}'.format(set(vertex)))
        if len(positions) < 2:
            raise ValueError('At least two vertices are needed')
        if not closed:
            assert len(tcb) >= 2
            if tcb.pop(0) != (0, 0, 0):
                raise ValueError(
                    'First vertex cannot have tension/continuity/bias (except for closed curves')
            if tcb.pop(-1) != (0, 0, 0):
                raise ValueError(
                    'Last vertex cannot have tension/continuity/bias (except for closed curves)')
        positions = _np.ascontiguousarray(positions, dtype='float32')
        assert positions.ndim == 2
        assert positions.shape[1] == 3
        positions_ptr = _ffi.from_buffer('float[]', positions)
        times = _np.ascontiguousarray(times, dtype='float32')
        assert times.ndim == 1
        times_ptr = _ffi.from_buffer('float[]', times)
        speeds = _np.ascontiguousarray(speeds, dtype='float32')
        assert speeds.ndim == 1
        speeds_ptr = _ffi.from_buffer('float[]', speeds)
        tcb = _np.ascontiguousarray(tcb, dtype='float32')
        tcb_ptr = _ffi.from_buffer('float[]', tcb)
        ptr = _ffi.gc(
            _lib.asdf_asdfspline(
                positions_ptr, len(positions),
                times_ptr, len(times),
                speeds_ptr, len(speeds),
                tcb_ptr, len(tcb),
                closed,
                ),
            _lib.asdf_asdfspline_free)
        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
        self._ptr = ptr

    def evaluate(self, t):
        t = _np.ascontiguousarray(t, dtype='float32')
        output = _np.empty(t.shape + (3,), dtype='float32')
        _lib.asdf_asdfspline_evaluate(
            self._ptr,
            _ffi.from_buffer('float[]', t),
            t.size,
            _ffi.from_buffer('float[]', output))
        return output

    @property
    def grid(self):
        grid_ptr = _ffi.new('float**')
        grid_len = _lib.asdf_asdfspline_grid(self._ptr, grid_ptr)
        if grid_len == 0:
            raise RuntimeError(_ffi.string(_lib.asdf_last_error()).decode())
        buffer = _ffi.buffer(
            grid_ptr[0], grid_len * _np.dtype('float32').itemsize)
        array = _np.frombuffer(buffer, dtype='float32')
        # NB: Writing to this array would be Undefined Behavior
        array.flags.writeable = False
        return array


class CubicCurve3:

    def __init__(self, ptr):
        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
        self._ptr = ptr

    def evaluate(self, t):
        t = _np.ascontiguousarray(t, dtype='float32')
        output = _np.empty(t.shape + (3,), dtype='float32')
        _lib.asdf_cubiccurve3_evaluate(
            self._ptr,
            _ffi.from_buffer('float[]', t),
            t.size,
            _ffi.from_buffer('float[]', output))
        return output

    @property
    def grid(self):
        grid_ptr = _ffi.new('float**')
        grid_len = _lib.asdf_cubiccurve3_grid(self._ptr, grid_ptr)
        if grid_len == 0:
            raise RuntimeError(_ffi.string(_lib.asdf_last_error()).decode())
        buffer = _ffi.buffer(
            grid_ptr[0], grid_len * _np.dtype('float32').itemsize)
        array = _np.frombuffer(buffer, dtype='float32')
        # NB: Writing to this array would be Undefined Behavior
        array.flags.writeable = False
        return array


class CentripetalKochanekBartelsSpline3(CubicCurve3):

    def __init__(self, positions, tcb=None, closed=False):
        positions = _np.ascontiguousarray(positions, dtype='float32')
        if positions.ndim != 2:
            raise ValueError(
                "Positions must be two-dimensional (list of coordinate pairs)")
        if positions.shape[1] != 3:
            raise ValueError("Positions must be a list of coordinate pairs")

        if tcb is None:
            tcb = 0, 0, 0
        tcb = _np.ascontiguousarray(tcb, dtype='float32')
        if tcb.ndim == 1:
            N = len(positions) if closed else len(positions) - 2
            tcb = _np.tile(tcb, (N, 1))
        if tcb.ndim != 2:
            raise ValueError(
                "TCB values must be two-dimensional (list of triples)")
        if tcb.shape[1] != 3:
            raise ValueError("TCB values must be a list of triples")

        ptr = _ffi.gc(
            _lib.asdf_centripetalkochanekbartelsspline3(
                _ffi.from_buffer('float[]', positions),
                len(positions),
                _ffi.from_buffer('float[]', tcb),
                len(tcb),
                closed,
            ),
            _lib.asdf_cubiccurve3_free)
        super().__init__(ptr)


class CubicCurve2:

    def __init__(self, ptr):
        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
        self._ptr = ptr

    def evaluate(self, t):
        t = _np.ascontiguousarray(t, dtype='float32')
        output = _np.empty(t.shape + (2,), dtype='float32')
        _lib.asdf_cubiccurve2_evaluate(
            self._ptr,
            _ffi.from_buffer('float[]', t),
            t.size,
            _ffi.from_buffer('float[]', output))
        return output

    @property
    def grid(self):
        grid_ptr = _ffi.new('float**')
        grid_len = _lib.asdf_cubiccurve2_grid(self._ptr, grid_ptr)
        if grid_len == 0:
            raise RuntimeError(_ffi.string(_lib.asdf_last_error()).decode())
        buffer = _ffi.buffer(
            grid_ptr[0], grid_len * _np.dtype('float32').itemsize)
        array = _np.frombuffer(buffer, dtype='float32')
        # NB: Writing to this array would be Undefined Behavior
        array.flags.writeable = False
        return array


class CentripetalKochanekBartelsSpline2(CubicCurve2):

    def __init__(self, positions, tcb=None, closed=False):
        positions = _np.ascontiguousarray(positions, dtype='float32')
        if positions.ndim != 2:
            raise ValueError(
                "Positions must be two-dimensional (list of coordinate pairs)")
        if positions.shape[1] != 2:
            raise ValueError("Positions must be a list of coordinate pairs")

        if tcb is None:
            tcb = 0, 0, 0
        tcb = _np.ascontiguousarray(tcb, dtype='float32')
        if tcb.ndim == 1:
            N = len(positions) if closed else len(positions) - 2
            tcb = _np.tile(tcb, (N, 1))
        if tcb.ndim != 2:
            raise ValueError(
                "TCB values must be two-dimensional (list of triples)")
        if tcb.shape[1] != 3:
            raise ValueError("TCB values must be a list of triples")

        ptr = _ffi.gc(
            _lib.asdf_centripetalkochanekbartelsspline2(
                _ffi.from_buffer('float[]', positions),
                len(positions),
                _ffi.from_buffer('float[]', tcb),
                len(tcb),
                closed,
            ),
            _lib.asdf_cubiccurve2_free)
        super().__init__(ptr)


class CubicCurve1:

    def __init__(self, ptr):
        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
        self._ptr = ptr

    def evaluate(self, t):
        t = _np.ascontiguousarray(t, dtype='float32')
        output = _np.empty_like(t)
        _lib.asdf_cubiccurve1_evaluate(
            self._ptr,
            _ffi.from_buffer('float[]', t),
            t.size,
            _ffi.from_buffer('float[]', output))
        return output

    @property
    def grid(self):
        grid_ptr = _ffi.new('float**')
        grid_len = _lib.asdf_cubiccurve1_grid(self._ptr, grid_ptr)
        if grid_len == 0:
            raise RuntimeError(_ffi.string(_lib.asdf_last_error()).decode())
        buffer = _ffi.buffer(
            grid_ptr[0], grid_len * _np.dtype('float32').itemsize)
        array = _np.frombuffer(buffer, dtype='float32')
        # NB: Writing to this array would be Undefined Behavior
        array.flags.writeable = False
        return array


class ShapePreservingCubicSpline(CubicCurve1):

    def __init__(self, values, slopes=None, grid=None, closed=False):
        values = _np.ascontiguousarray(values, dtype='float32')
        if values.ndim != 1:
            raise ValueError("Values must be one-dimensional")
        if grid is None:
            grid = _np.arange(len(values) + closed, dtype='float32')
        grid = _np.ascontiguousarray(grid, dtype='float32')
        if grid.ndim != 1:
            raise ValueError("Grid must be one-dimensional")
        values_ptr = _ffi.from_buffer('float[]', values)
        grid_ptr = _ffi.from_buffer('float[]', grid)
        if slopes is None:
            ptr = _ffi.gc(
                _lib.asdf_shapepreservingcubicspline(
                    values_ptr, len(values),
                    grid_ptr, len(grid),
                    closed,
                ),
                _lib.asdf_cubiccurve1_free)
        else:
            slopes = _np.ascontiguousarray(slopes, dtype='float32')
            if slopes.ndim != 1:
                raise ValueError("Slopes must be one-dimensional")
            slopes_ptr = _ffi.from_buffer('float[]', slopes)
            ptr = _ffi.gc(
                _lib.asdf_shapepreservingcubicspline_with_slopes(
                    values_ptr, len(values),
                    slopes_ptr, len(slopes),
                    grid_ptr, len(grid),
                    closed,
                ),
                _lib.asdf_cubiccurve1_free)
        super().__init__(ptr)


class MonotoneCubicSpline(CubicCurve1):

    def __init__(self, values, slopes=None, grid=None):
        # TODO: code re-use?
        values = _np.ascontiguousarray(values, dtype='float32')
        if values.ndim != 1:
            raise ValueError("Values must be one-dimensional")
        if grid is None:
            grid = _np.arange(len(values), dtype='float32')
        grid = _np.ascontiguousarray(grid, dtype='float32')
        if grid.ndim != 1:
            raise ValueError("Grid must be one-dimensional")
        values_ptr = _ffi.from_buffer('float[]', values)
        grid_ptr = _ffi.from_buffer('float[]', grid)
        if slopes is None:
            ptr = _ffi.gc(
                _lib.asdf_monotonecubic(
                    values_ptr, len(values),
                    grid_ptr, len(grid),
                ),
                _lib.asdf_monotonecubic_free)
        else:
            # TODO: code re-use?
            slopes = _np.ascontiguousarray(slopes, dtype='float32')
            if slopes.ndim != 1:
                raise ValueError("Slopes must be one-dimensional")
            slopes_ptr = _ffi.from_buffer('float[]', slopes)
            ptr = _ffi.gc(
                _lib.asdf_monotonecubic_with_slopes(
                    values_ptr, len(values),
                    slopes_ptr, len(slopes),
                    grid_ptr, len(grid),
                ),
                _lib.asdf_monotonecubic_free)

        if ptr == _ffi.NULL:
            raise ValueError(_ffi.string(_lib.asdf_last_error()).decode())
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
