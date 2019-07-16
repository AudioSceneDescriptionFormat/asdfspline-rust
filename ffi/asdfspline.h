#ifndef ASDFSPLINE_H
#define ASDFSPLINE_H

/* Generated with cbindgen:0.9.0 */

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct AsdfSpline_f32__Vec3 AsdfSpline_f32__Vec3;

/**
 * ... monotonically *increasing* ...
 */
typedef struct MonotoneCubicSpline_f32 MonotoneCubicSpline_f32;

typedef struct PiecewiseCubicCurve_f32__Vec2 PiecewiseCubicCurve_f32__Vec2;

typedef struct PiecewiseCubicCurve_f32__Vec3 PiecewiseCubicCurve_f32__Vec3;

typedef struct PiecewiseCubicCurve_f32__f32 PiecewiseCubicCurve_f32__f32;

typedef AsdfSpline_f32__Vec3 AsdfSpline;

typedef PiecewiseCubicCurve_f32__Vec2 AsdfCubicCurve2;

typedef PiecewiseCubicCurve_f32__Vec3 AsdfCubicCurve3;

typedef PiecewiseCubicCurve_f32__f32 AsdfCubicCurve1;

typedef MonotoneCubicSpline_f32 AsdfMonotoneCubic;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

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

/**
 * The error message will be freed if another error occurs. It is the caller's
 * responsibility to make sure they're no longer using the string before
 * calling any other function which may fail.
 */
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

#ifdef __cplusplus
} // extern "C"
#endif // __cplusplus

#endif /* ASDFSPLINE_H */
