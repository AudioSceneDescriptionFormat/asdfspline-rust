from cffi import FFI

from subprocess import run

run(
    ['cargo', 'build', '--all', '--release'],
    cwd='..',
    check=True,
)

ffibuilder = FFI()
ffibuilder.cdef("""
/* Opaque structs: */
typedef struct AsdfSpline AsdfSpline;
typedef struct AsdfCubicCurve3 AsdfCubicCurve3;
typedef struct AsdfCubicCurve2 AsdfCubicCurve2;
typedef struct AsdfCubicCurve1 AsdfCubicCurve1;
typedef struct AsdfMonotoneCubic AsdfMonotoneCubic;

/* Copied from ffi/asdfspline.h: */
AsdfSpline *asdf_asdfspline(const float *positions,
                            size_t positions_count,
                            const float *times,
                            size_t times_count,
                            const float *speeds,
                            size_t speeds_count,
                            const float *tcb,
                            size_t tcb_count,
                            bool closed);
void asdf_asdfspline_evaluate(AsdfSpline *ptr, const float *times, size_t count, float *output);
void asdf_asdfspline_free(AsdfSpline *ptr);
size_t asdf_asdfspline_grid(AsdfSpline *ptr, const float **output);
AsdfCubicCurve2 *asdf_centripetalkochanekbartelsspline2(const float *vertices,
                                                        size_t vertices_count,
                                                        const float *tcb,
                                                        size_t tcb_count,
                                                        bool closed);
AsdfCubicCurve3 *asdf_centripetalkochanekbartelsspline3(const float *vertices,
                                                        size_t vertices_count,
                                                        const float *tcb,
                                                        size_t tcb_count,
                                                        bool closed);
void asdf_cubiccurve1_evaluate(AsdfCubicCurve1 *ptr,
                               const float *times,
                               size_t count,
                               float *output);
void asdf_cubiccurve1_free(AsdfCubicCurve1 *ptr);
size_t asdf_cubiccurve1_grid(AsdfCubicCurve1 *ptr, const float **output);
void asdf_cubiccurve2_evaluate(AsdfCubicCurve2 *ptr,
                               const float *times,
                               size_t count,
                               float *output);
void asdf_cubiccurve2_free(AsdfCubicCurve2 *ptr);
size_t asdf_cubiccurve2_grid(AsdfCubicCurve2 *ptr, const float **output);
void asdf_cubiccurve3_evaluate(AsdfCubicCurve3 *ptr,
                               const float *times,
                               size_t count,
                               float *output);
void asdf_cubiccurve3_free(AsdfCubicCurve3 *ptr);
size_t asdf_cubiccurve3_grid(AsdfCubicCurve3 *ptr, const float **output);
const char *asdf_last_error(void);
AsdfMonotoneCubic *asdf_monotonecubic(const float *values,
                                      size_t values_count,
                                      const float *grid,
                                      size_t grid_count);
void asdf_monotonecubic_free(AsdfMonotoneCubic *ptr);
void asdf_monotonecubic_get_time(AsdfMonotoneCubic *ptr,
                                 const float *values,
                                 size_t count,
                                 float *output);
const AsdfCubicCurve1 *asdf_monotonecubic_inner(AsdfMonotoneCubic *ptr);
AsdfMonotoneCubic *asdf_monotonecubic_with_slopes(const float *values,
                                                  size_t values_count,
                                                  const float *slopes,
                                                  size_t slopes_count,
                                                  const float *grid,
                                                  size_t grid_count);
AsdfCubicCurve1 *asdf_shapepreservingcubicspline(const float *values,
                                                 size_t values_count,
                                                 const float *grid,
                                                 size_t grid_count,
                                                 bool closed);
AsdfCubicCurve1 *asdf_shapepreservingcubicspline_with_slopes(const float *values,
                                                             size_t values_count,
                                                             const float *slopes,
                                                             size_t slopes_count,
                                                             const float *grid,
                                                             size_t grid_count,
                                                             bool closed);
""")
ffibuilder.set_source(
    '_asdfspline',
    r"""
    #include "asdfspline.h"
    """,
    include_dirs=['../ffi'],
    libraries=['asdfspline_ffi'],
    library_dirs=['../target/release'],
    undef_macros=['NDEBUG'],
)

if __name__ == '__main__':
    ffibuilder.compile(verbose=True)
