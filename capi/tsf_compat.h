/* tsf_compat.h -- C API for reading Bruker timsTOF TSF data.
 *
 * This header declares the same function signatures as the vendor SDK's TSF
 * "mini API", so that software written against that API can load this library
 * instead on platforms the vendor does not ship a build for -- macOS in
 * particular. It is an independent declaration of the interface, written from
 * the documented signatures; no vendor header or source is reproduced here.
 *
 * Part of pytsf. MIT licensed. See LICENSE and PROVENANCE.md.
 *
 * Differences from the vendor library, which callers should know about:
 *
 *   - tsf_has_recalibrated_state() always returns 0: post-hoc recalibrated
 *     states are not supported. Raw acquisition-time calibration is used, and
 *     tsf_open() fails if a recalibrated state is requested.
 *   - tsf_set_num_threads() is accepted and ignored; this implementation is
 *     single-threaded per call, and callers may use separate handles from
 *     separate threads.
 *   - Passing length <= 0 to a read function is a pure size query: nothing is
 *     written, whatever the array pointers are, and the required count is
 *     returned.
 *   - m/z is obtained by fitting the published three-point time-of-flight
 *     equation -- t = ml2 + sqrt(1e12/ml1 * m) + ml3*m, per OpenMS
 *     TOFCalibration after Gobom et al., Anal. Chem. 2002, 74(15), 3915-3923 --
 *     to the reference masses and measured flight times the instrument
 *     recorded in the file. That is an external calibration: accurate across
 *     the calibrant's mass range and degrading outside it. Expect agreement
 *     with the vendor library of order 1 ppm, not exactness. tsf_open() fails
 *     if the file records no calibration points.
 *   - Only the TSF format is implemented here. TDF (ion mobility) is read by the
 *     Python package (pytsf.TdfData); a tims_* C ABI is not yet provided.
 */

#ifndef PYTSF_TSF_COMPAT_H
#define PYTSF_TSF_COMPAT_H

#include <stdint.h>

#ifdef _WIN32
#  define TSF_API __declspec(dllexport)
#else
#  define TSF_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Open an analysis directory (UTF-8). Returns a non-zero handle on success,
 * 0 on failure (use tsf_get_last_error_string). use_recalibrated_state must
 * be 0; a non-zero value fails rather than silently returning raw values. */
TSF_API uint64_t tsf_open(const char *analysis_directory_name,
                          uint32_t use_recalibrated_state);

/* Close a handle. Passing 0 is safe and does nothing. */
TSF_API void tsf_close(uint64_t handle);

/* Last error for the calling thread. Returns the full length including the
 * terminating zero; if that exceeds len the message was truncated. */
TSF_API uint32_t tsf_get_last_error_string(char *buf, uint32_t len);

/* Always 0 in this implementation (see header comment). */
TSF_API uint32_t tsf_has_recalibrated_state(uint64_t handle);

/* Read a centroided ("line") spectrum. Returns -1 on error, else the number
 * of entries the output arrays need; if that exceeds `length` the result is
 * incomplete and the call should be repeated with a larger buffer. */
TSF_API int32_t tsf_read_line_spectrum_v2(uint64_t handle, int64_t spectrum_id,
                                          double *index_array,
                                          float *intensity_array,
                                          int32_t length);

/* As above, additionally returning peak widths. Fails if the analysis stores
 * no peak widths. */
TSF_API int32_t tsf_read_line_spectrum_with_width_v2(uint64_t handle,
                                                     int64_t spectrum_id,
                                                     double *index_array,
                                                     float *intensity_array,
                                                     float *width_array,
                                                     int32_t length);

/* Read a profile spectrum into a dense array; position in the array is the
 * index. Same return convention as above. */
TSF_API int32_t tsf_read_profile_spectrum_v2(uint64_t handle,
                                             int64_t spectrum_id,
                                             uint32_t *profile_array,
                                             int32_t length);

/* Conversion functions. Return 1 on success, 0 on failure. */
TSF_API uint32_t tsf_index_to_mz(uint64_t handle, int64_t frame_id,
                                 const double *in, double *out, uint32_t cnt);
TSF_API uint32_t tsf_mz_to_index(uint64_t handle, int64_t frame_id,
                                 const double *in, double *out, uint32_t cnt);

/* Accepted and ignored. */
TSF_API void tsf_set_num_threads(uint32_t n);

/* --- extensions beyond the vendor API ------------------------------------ */

/* Library version string, e.g. "pytsf-capi 1.0.0". Never NULL. */
TSF_API const char *tsf_compat_version(void);

#ifdef __cplusplus
}
#endif

#endif /* PYTSF_TSF_COMPAT_H */
