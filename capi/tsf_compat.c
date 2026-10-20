/* tsf_compat.c -- reading Bruker timsTOF TSF data, in C, with no vendor code.
 *
 * Part of pytsf. MIT licensed. See LICENSE and PROVENANCE.md for where the
 * understanding of the format comes from and how this is validated.
 *
 * Depends only on SQLite and zstd. Mirrors the Python implementation in
 * pytsf/bin.py and pytsf/calibration_fitted.py; the two are validated against
 * each other and against the vendor library.
 */

/* Must precede every system header so that off_t / fseeko are 64-bit on
 * 32-bit-long platforms. */
#if !defined(_WIN32)
#  define _FILE_OFFSET_BITS 64
#  define _DEFAULT_SOURCE 1
#endif

#include "tsf_compat.h"

#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sqlite3.h>
#include <zstd.h>

/* 64-bit file positioning: `long` is 32 bits on Windows and imaging .tsf_bin
 * files routinely exceed 2 GiB. (_FILE_OFFSET_BITS is set at the top of this
 * file, before any system header.) */
#if defined(_WIN32)
#  define tsf_fseek(f, off) _fseeki64((f), (long long)(off), SEEK_SET)
#else
#  define tsf_fseek(f, off) fseeko((f), (off_t)(off), SEEK_SET)
#endif

#define TSF_COMPAT_VERSION "pytsf-capi 1.0.0"

/* Largest sub-block we will attempt to read (1 GiB), matching the vendor
 * library's limit. */
#define MAX_SUBBLOCK (0x40000000u)

/* ------------------------------------------------------------------ errors */

#if defined(_MSC_VER)
#  define TSF_THREAD __declspec(thread)
#else
#  define TSF_THREAD __thread
#endif

static TSF_THREAD char g_error[1024];

static void set_error(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(g_error, sizeof(g_error), fmt, ap);
    va_end(ap);
}

uint32_t tsf_get_last_error_string(char *buf, uint32_t len) {
    uint32_t need = (uint32_t)strlen(g_error) + 1u;
    if (buf && len) {
        uint32_t n = need < len ? need : len;
        memcpy(buf, g_error, n);
        buf[n - 1] = '\0';
    }
    return need;
}

/* ------------------------------------------------- fitted calibration ----
 *
 * The published three-point time-of-flight equation
 *
 *     t = ml2 + sqrt(1e12/ml1 * m) + ml3 * m
 *
 * (OpenMS TOFCalibration; after Gobom et al., Anal. Chem. 2002, 74(15),
 * 3915-3923), fitted by Gauss-Newton to the reference masses and measured
 * flight times the instrument recorded in the documented CalibrationInfo
 * table. It rests only on published and vendor-documented material; see
 * PROVENANCE.md.
 */

#define MAX_REFPTS 64
#define NPARAM 3

typedef struct {
    int    valid;
    double p[NPARAM];        /* ml1, ml2, ml3 */
    double timebase, delay;
    double mz_min, mz_max;   /* calibrant range; outside it we extrapolate */
    int    npts;
    double fit_rms_ppm;
} FittedCalib;

static double fitted_forward(const double *p, double m) {
    double am = fabs(m);
    return p[1] + sqrt(1e12 / p[0] * am) + p[2] * m;
}

/* Flight times at or before the fitted time origin have no physical solution
 * -- no ion arrives before the start of flight -- so return NaN rather than a
 * fabricated value. */
static double fitted_mz_from_time(const FittedCalib *fc, double t) {
    double d = t - fc->p[1];
    double m;
    int i;
    if (!(d > 0.0)) return NAN;
    m = d * d / (1e12 / fc->p[0]);
    for (i = 0; i < 80; i++) {
        double h  = fabs(m) * 1e-7; double f, fp, mn;
        if (h < 1e-10) h = 1e-10;
        f  = fitted_forward(fc->p, m) - t;
        fp = (fitted_forward(fc->p, m + h) - fitted_forward(fc->p, m - h)) / (2.0 * h);
        if (fp == 0.0) break;
        mn = m - f / fp;
        if (!(mn > 0.0)) mn = m * 0.5;
        if (fabs(mn - m) <= 1e-13 * fabs(m)) { m = mn; break; }
        m = mn;
    }
    return m;
}

/* Solve the small normal equations for the Gauss-Newton step. */
static int solve_normal(double JtJ[NPARAM][NPARAM], double Jtr[NPARAM],
                        double out[NPARAM]) {
    int i, j, k, piv;
    double A[NPARAM][NPARAM + 1];
    for (i = 0; i < NPARAM; i++) {
        for (j = 0; j < NPARAM; j++) A[i][j] = JtJ[i][j];
        A[i][NPARAM] = Jtr[i];
    }
    for (k = 0; k < NPARAM; k++) {
        double best = fabs(A[k][k]); piv = k;
        for (i = k + 1; i < NPARAM; i++)
            if (fabs(A[i][k]) > best) { best = fabs(A[i][k]); piv = i; }
        if (best < 1e-300) return 0;
        if (piv != k) for (j = 0; j <= NPARAM; j++) {
            double tmp = A[k][j]; A[k][j] = A[piv][j]; A[piv][j] = tmp;
        }
        for (i = k + 1; i < NPARAM; i++) {
            double f = A[i][k] / A[k][k];
            for (j = k; j <= NPARAM; j++) A[i][j] -= f * A[k][j];
        }
    }
    for (i = NPARAM - 1; i >= 0; i--) {
        double sum = A[i][NPARAM];
        for (j = i + 1; j < NPARAM; j++) sum -= A[i][j] * out[j];
        out[i] = sum / A[i][i];
    }
    return 1;
}

static int fitted_fit(FittedCalib *fc, const double *m, const double *t, int n) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0, slope, intercept, det;
    int i, it;

    /* seed: t is affine in sqrt(m) */
    for (i = 0; i < n; i++) {
        double x = sqrt(m[i]);
        sx += x; sy += t[i]; sxx += x * x; sxy += x * t[i];
    }
    det = n * sxx - sx * sx;
    if (fabs(det) < 1e-300) { set_error("Degenerate calibration points."); return 0; }
    slope     = (n * sxy - sx * sy) / det;
    intercept = (sy - slope * sx) / n;
    if (slope <= 0.0) { set_error("Calibration points are not physically ordered."); return 0; }
    fc->p[0] = 1e12 / (slope * slope);
    fc->p[1] = intercept;
    fc->p[2] = 0.0;

    for (it = 0; it < 100; it++) {
        double JtJ[NPARAM][NPARAM], Jtr[NPARAM], step[NPARAM];
        int a, b, converged = 1;
        for (a = 0; a < NPARAM; a++) {
            Jtr[a] = 0.0;
            for (b = 0; b < NPARAM; b++) JtJ[a][b] = 0.0;
        }
        for (i = 0; i < n; i++) {
            double J[NPARAM], r = fitted_forward(fc->p, m[i]) - t[i];
            for (a = 0; a < NPARAM; a++) {
                double save = fc->p[a];
                double h = fabs(save) * 1e-7; double f1;
                if (h < 1e-10) h = 1e-10;
                fc->p[a] = save + h;
                f1 = fitted_forward(fc->p, m[i]);
                fc->p[a] = save;
                J[a] = (f1 - fitted_forward(fc->p, m[i])) / h;
            }
            for (a = 0; a < NPARAM; a++) {
                Jtr[a] -= J[a] * r;
                for (b = 0; b < NPARAM; b++) JtJ[a][b] += J[a] * J[b];
            }
        }
        if (!solve_normal(JtJ, Jtr, step)) break;
        for (a = 0; a < NPARAM; a++) {
            if (fabs(step[a]) > 1e-14 * (fabs(fc->p[a]) > 1.0 ? fabs(fc->p[a]) : 1.0))
                converged = 0;
            fc->p[a] += step[a];
        }
        if (!(fc->p[0] > 0.0) || !isfinite(fc->p[0])) {
            set_error("Calibration fit diverged."); return 0;
        }
        if (converged) break;
    }

    {
        double acc = 0.0;
        fc->mz_min = fc->mz_max = m[0];
        for (i = 0; i < n; i++) {
            double got = fitted_mz_from_time(fc, t[i]);
            double rel = (got - m[i]) / m[i];
            acc += rel * rel;
            if (m[i] < fc->mz_min) fc->mz_min = m[i];
            if (m[i] > fc->mz_max) fc->mz_max = m[i];
        }
        fc->fit_rms_ppm = sqrt(acc / n) * 1e6;
    }
    fc->npts = n;
    fc->valid = 1;
    return 1;
}

/* ---------------------------------------------------------------- instance */

typedef struct {
    sqlite3 *db;
    FILE    *bin;
    FittedCalib fitted;
    int      has_line, has_profile, has_width;
    /* scratch buffers */
    unsigned char *comp;   size_t comp_cap;
    unsigned char *raw;    size_t raw_cap;
} Instance;

static int ensure_cap(unsigned char **buf, size_t *cap, size_t need) {
    if (*cap >= need) return 1;
    {
        unsigned char *p = (unsigned char *)realloc(*buf, need);
        if (!p) return 0;
        *buf = p; *cap = need;
    }
    return 1;
}

static int meta_int(sqlite3 *db, const char *key, long long *out) {
    sqlite3_stmt *st = NULL;
    int rc, got = 0;
    if (sqlite3_prepare_v2(db, "SELECT Value FROM GlobalMetadata WHERE Key=?",
                           -1, &st, NULL) != SQLITE_OK) return 0;
    sqlite3_bind_text(st, 1, key, -1, SQLITE_STATIC);
    rc = sqlite3_step(st);
    if (rc == SQLITE_ROW) {
        const unsigned char *v = sqlite3_column_text(st, 0);
        if (v) { *out = atoll((const char *)v); got = 1; }
    }
    sqlite3_finalize(st);
    return got;
}

/* Read a BLOB of little-endian doubles from CalibrationInfo. Returns count,
 * or -1 if absent/unusable. */
static int calinfo_doubles(sqlite3 *db, const char *key, double *out, int cap) {
    sqlite3_stmt *st = NULL;
    int n = -1;
    if (sqlite3_prepare_v2(db, "SELECT Value FROM CalibrationInfo WHERE KeyName=?",
                           -1, &st, NULL) != SQLITE_OK) return -1;
    sqlite3_bind_text(st, 1, key, -1, SQLITE_STATIC);
    if (sqlite3_step(st) == SQLITE_ROW
        && sqlite3_column_type(st, 0) == SQLITE_BLOB) {
        int bytes = sqlite3_column_bytes(st, 0);
        const void *blob = sqlite3_column_blob(st, 0);
        if (blob && bytes >= 8 && bytes % 8 == 0) {
            n = bytes / 8;
            if (n > cap) {
                n = -2;   /* too many points: caller reports it */
            } else {
                memcpy(out, blob, (size_t)n * 8);
            }
        }
    }
    sqlite3_finalize(st);
    return n;
}

static int load_fitted(Instance *inst) {
    double mz[MAX_REFPTS], tof[MAX_REFPTS], m2[MAX_REFPTS], t2[MAX_REFPTS];
    int nm, nt, n, i, k = 0;
    sqlite3_stmt *st = NULL;

    nm = calinfo_doubles(inst->db, "ReferencePeakMasses", mz, MAX_REFPTS);
    nt = calinfo_doubles(inst->db, "MeasuredTimesOfFlight", tof, MAX_REFPTS);
    if (nm == -2 || nt == -2) {
        set_error("More than %d calibration points; raise MAX_REFPTS.",
                  MAX_REFPTS);
        return 0;
    }
    if (nm < 0 || nt < 0) {
        set_error("This analysis records no ReferencePeakMasses / "
                  "MeasuredTimesOfFlight in CalibrationInfo, so the calibration "
                  "cannot be fitted from the file.");
        return 0;
    }
    n = nm < nt ? nm : nt;
    for (i = 0; i < n; i++) {
        if (isfinite(mz[i]) && isfinite(tof[i]) && mz[i] > 0.0 && tof[i] > 0.0) {
            m2[k] = mz[i]; t2[k] = tof[i]; k++;
        }
    }
    if (k < NPARAM) {
        set_error("Only %d usable calibration points; need at least %d.",
                  k, NPARAM);
        return 0;
    }

    /* One fit covers the whole analysis, so it must have one digitiser
     * timing. Refuse a file whose spectra reference different calibration
     * rows rather than silently applying the first row's timing to all. */
    if (sqlite3_prepare_v2(inst->db,
            "SELECT COUNT(DISTINCT MzCalibration) FROM Frames",
            -1, &st, NULL) == SQLITE_OK && sqlite3_step(st) == SQLITE_ROW) {
        int ncal = sqlite3_column_int(st, 0);
        sqlite3_finalize(st); st = NULL;
        if (ncal > 1) {
            set_error("Spectra reference %d different MzCalibration rows; "
                      "this library fits one calibration per analysis and "
                      "does not support that.", ncal);
            return 0;
        }
    } else if (st) { sqlite3_finalize(st); st = NULL; }

    if (sqlite3_prepare_v2(inst->db,
            "SELECT DigitizerTimebase, DigitizerDelay FROM MzCalibration "
            "WHERE Id=(SELECT MzCalibration FROM Frames ORDER BY Id LIMIT 1)",
            -1, &st, NULL) != SQLITE_OK
        || sqlite3_step(st) != SQLITE_ROW) {
        if (st) sqlite3_finalize(st);
        set_error("No MzCalibration row found.");
        return 0;
    }
    inst->fitted.timebase = sqlite3_column_double(st, 0);
    inst->fitted.delay    = sqlite3_column_double(st, 1);
    sqlite3_finalize(st);
    if (inst->fitted.timebase == 0.0) {
        set_error("DigitizerTimebase must be non-zero."); return 0;
    }
    return fitted_fit(&inst->fitted, m2, t2, k);
}

uint64_t tsf_open(const char *dir, uint32_t use_recalibrated_state) {
    Instance *inst = NULL;
    char path[4096];
    long long v = 0;

    g_error[0] = '\0';
    if (!dir) { set_error("analysis_directory_name is NULL."); return 0; }
    if (use_recalibrated_state) {
        set_error("use_recalibrated_state is not supported by this "
                  "implementation; pass 0 to read the raw calibration.");
        return 0;
    }

    inst = (Instance *)calloc(1, sizeof(Instance));
    if (!inst) { set_error("Out of memory."); return 0; }

    snprintf(path, sizeof(path), "%s/analysis.tsf", dir);
    if (sqlite3_open_v2(path, &inst->db, SQLITE_OPEN_READONLY, NULL) != SQLITE_OK) {
        set_error("Cannot open %s: %s", path,
                  inst->db ? sqlite3_errmsg(inst->db) : "unknown error");
        goto fail;
    }

    {
        sqlite3_stmt *st = NULL;
        if (sqlite3_prepare_v2(inst->db,
                "SELECT Value FROM GlobalMetadata WHERE Key='SchemaType'",
                -1, &st, NULL) == SQLITE_OK) {
            if (sqlite3_step(st) == SQLITE_ROW) {
                const unsigned char *t = sqlite3_column_text(st, 0);
                if (t && strcmp((const char *)t, "TSF") != 0
                      && strcmp((const char *)t, "tsf") != 0) {
                    set_error("Not a TSF analysis: SchemaType=%s. This "
                              "library covers TSF only; for TDF use the Python package.", t);
                    sqlite3_finalize(st);
                    goto fail;
                }
            }
            sqlite3_finalize(st);
        }
    }

    if (!meta_int(inst->db, "TimsCompressionType", &v)) {
        set_error("GlobalMetadata is missing TimsCompressionType.");
        goto fail;
    }
    if (v != 3) {
        set_error("TimsCompressionType %lld not supported.", v);
        goto fail;
    }

    inst->has_line    = meta_int(inst->db, "HasLineSpectra", &v) ? (v != 0) : 0;
    inst->has_profile = meta_int(inst->db, "HasProfileSpectra", &v) ? (v != 0) : 0;
    inst->has_width   = meta_int(inst->db, "HasLineSpectraPeakWidth", &v) ? (v != 0) : 0;

    snprintf(path, sizeof(path), "%s/analysis.tsf_bin", dir);
    inst->bin = fopen(path, "rb");
    if (!inst->bin) { set_error("Cannot open %s", path); goto fail; }

    if (!load_fitted(inst)) goto fail;

    return (uint64_t)(uintptr_t)inst;

fail:
    if (inst->db) sqlite3_close(inst->db);
    if (inst->bin) fclose(inst->bin);
    free(inst);
    return 0;
}

void tsf_close(uint64_t handle) {
    Instance *inst = (Instance *)(uintptr_t)handle;
    if (!inst) return;
    if (inst->db) sqlite3_close(inst->db);
    if (inst->bin) fclose(inst->bin);
    free(inst->comp);
    free(inst->raw);
    free(inst);
}

uint32_t tsf_has_recalibrated_state(uint64_t handle) { (void)handle; return 0; }
void tsf_set_num_threads(uint32_t n) { (void)n; }
const char *tsf_compat_version(void) { return TSF_COMPAT_VERSION; }

/* ------------------------------------------------------- frame / calib I/O */

/* Look up TimsId / NumPeaks for one spectrum. */
static int frame_row(Instance *inst, int64_t fid, long long *tims_id,
                     long long *num_peaks, int *tims_null, int *peaks_null) {
    sqlite3_stmt *st = NULL;
    int rc, ok = 0;
    if (sqlite3_prepare_v2(inst->db,
            "SELECT TimsId, NumPeaks FROM Frames WHERE Id=?", -1, &st, NULL)
        != SQLITE_OK) {
        set_error("SQL error: %s", sqlite3_errmsg(inst->db));
        return 0;
    }
    sqlite3_bind_int64(st, 1, fid);
    rc = sqlite3_step(st);
    if (rc == SQLITE_ROW) {
        *tims_null  = sqlite3_column_type(st, 0) == SQLITE_NULL;
        *peaks_null = sqlite3_column_type(st, 1) == SQLITE_NULL;
        *tims_id    = *tims_null  ? -1 : sqlite3_column_int64(st, 0);
        *num_peaks  = *peaks_null ? -1 : sqlite3_column_int64(st, 1);
        ok = 1;
    } else {
        set_error("No such spectrum id: %lld", (long long)fid);
    }
    sqlite3_finalize(st);
    return ok;
}


/* --------------------------------------------------------------- blob I/O */

/* Read and decompress one sub-block. which: 0 = line, 1 = profile.
 * On success points raw/rawlen at inst->raw. Returns 1 ok, 0 error,
 * -1 "no data" (empty, not an error). */
static int read_subblock(Instance *inst, long long off, int which,
                         unsigned char **raw, size_t *rawlen) {
    unsigned char hdr[8];
    uint32_t total, line_size;
    long long start;
    size_t size;
    unsigned long long dsize;

    if (off < 0) return -1;
    if (tsf_fseek(inst->bin, off) != 0) {
        set_error("Cannot seek to binary header at offset %lld", off); return 0;
    }
    if (fread(hdr, 1, 8, inst->bin) != 8) {
        set_error("Error reading binary header at offset %lld", off); return 0;
    }
    total     = (uint32_t)hdr[0] | ((uint32_t)hdr[1] << 8)
              | ((uint32_t)hdr[2] << 16) | ((uint32_t)hdr[3] << 24);
    line_size = (uint32_t)hdr[4] | ((uint32_t)hdr[5] << 8)
              | ((uint32_t)hdr[6] << 16) | ((uint32_t)hdr[7] << 24);
    if (total < 8u || line_size > total - 8u) {
        set_error("Corrupt block header at offset %lld", off); return 0;
    }

    if (which == 0) { start = off + 8;             size = line_size; }
    else            { start = off + 8 + line_size; size = total - 8u - line_size; }
    if (size == 0) return -1;
    if (size > MAX_SUBBLOCK) {
        set_error("Compressed spectrum size too large (%llu).",
                  (unsigned long long)size);
        return 0;
    }

    if (!ensure_cap(&inst->comp, &inst->comp_cap, size)) {
        set_error("Out of memory."); return 0;
    }
    if (tsf_fseek(inst->bin, start) != 0 ||
        fread(inst->comp, 1, size, inst->bin) != size) {
        set_error("Error reading frame data at offset %lld", (long long)start);
        return 0;
    }

    dsize = ZSTD_getFrameContentSize(inst->comp, size);
    if (dsize == ZSTD_CONTENTSIZE_ERROR) {
        set_error("Error decompressing frame data: not a zstd frame."); return 0;
    }
    if (dsize == ZSTD_CONTENTSIZE_UNKNOWN) {
        set_error("Error decompressing frame data: unknown content size.");
        return 0;
    }
    if (dsize == 0) return -1;
    if (!ensure_cap(&inst->raw, &inst->raw_cap, (size_t)dsize)) {
        set_error("Out of memory."); return 0;
    }
    {
        size_t got = ZSTD_decompress(inst->raw, (size_t)dsize, inst->comp, size);
        if (ZSTD_isError(got)) {
            set_error("Error decompressing frame data: %s",
                      ZSTD_getErrorName(got));
            return 0;
        }
        *raw = inst->raw; *rawlen = got;
    }
    return 1;
}

static void le_f64(const unsigned char *p, double *out, size_t n) {
    memcpy(out, p, n * 8);   /* x86-64 and arm64 are both little-endian */
}
static void le_f32(const unsigned char *p, float *out, size_t n) {
    memcpy(out, p, n * 4);
}

static int32_t read_line(uint64_t handle, int64_t fid, double *idx, float *inten,
                         float *width, int32_t length, int need_width) {
    Instance *inst = (Instance *)(uintptr_t)handle;
    long long tims_id = -1, num_peaks = -1;
    int tnull = 0, pnull = 0, rc;
    unsigned char *raw = NULL; size_t rawlen = 0;
    size_t n, need, have;

    g_error[0] = '\0';
    if (!inst) { set_error("Invalid handle."); return -1; }
    if (!inst->has_line) {
        set_error("This analysis contains no line spectra "
                  "(GlobalMetadata.HasLineSpectra == 0).");
        return -1;
    }
    if (need_width && !inst->has_width) {
        set_error("This analysis contains no peak widths "
                  "(GlobalMetadata.HasLineSpectraPeakWidth == 0).");
        return -1;
    }
    if (!frame_row(inst, fid, &tims_id, &num_peaks, &tnull, &pnull)) return -1;
    if (tnull || pnull || num_peaks <= 0) return 0;

    rc = read_subblock(inst, tims_id, 0, &raw, &rawlen);
    if (rc == 0) return -1;
    if (rc < 0) return 0;

    n = (size_t)num_peaks;
    /* Infer the width plane from the payload length rather than trusting the
     * metadata flag: a mismatch would otherwise read out of bounds. */
    have = rawlen >= n * 16 ? 16 : 12;
    if (need_width && have < 16) {
        set_error("Spectrum %lld stores no peak widths.", (long long)fid);
        return -1;
    }
    need = n * have;
    if (rawlen < need) {
        set_error("Line sub-block too short: %llu bytes, need %llu for %llu peaks",
                  (unsigned long long)rawlen, (unsigned long long)need,
                  (unsigned long long)n);
        return -1;
    }
    /* Copy at most `length` entries. length <= 0 is a pure size query and
     * must write nothing, whatever the pointers are: callers commonly pass a
     * zero-length buffer to discover the size. The return value is always the
     * full count the caller needs. */
    {
        size_t m = (length > 0 && (size_t)length < n) ? (size_t)length
                 : (length > 0 ? n : 0);
        if (m) {
            if (idx)   le_f64(raw, idx, m);
            if (inten) le_f32(raw + n * 8, inten, m);
            if (width && have == 16) le_f32(raw + n * 12, width, m);
        }
    }
    return (int32_t)n;
}

int32_t tsf_read_line_spectrum_v2(uint64_t h, int64_t fid, double *idx,
                                  float *inten, int32_t length) {
    return read_line(h, fid, idx, inten, NULL, length, 0);
}

int32_t tsf_read_line_spectrum_with_width_v2(uint64_t h, int64_t fid, double *idx,
                                             float *inten, float *width,
                                             int32_t length) {
    return read_line(h, fid, idx, inten, width, length, 1);
}

int32_t tsf_read_profile_spectrum_v2(uint64_t handle, int64_t fid,
                                     uint32_t *profile, int32_t length) {
    Instance *inst = (Instance *)(uintptr_t)handle;
    long long tims_id = -1, num_peaks = -1;
    int tnull = 0, pnull = 0, rc;
    unsigned char *raw = NULL; size_t rawlen = 0;
    size_t n, i, m;

    g_error[0] = '\0';
    if (!inst) { set_error("Invalid handle."); return -1; }
    if (!inst->has_profile) {
        set_error("This analysis contains no profile spectra "
                  "(GlobalMetadata.HasProfileSpectra == 0).");
        return -1;
    }
    if (!frame_row(inst, fid, &tims_id, &num_peaks, &tnull, &pnull)) return -1;
    if (tnull) return 0;

    rc = read_subblock(inst, tims_id, 1, &raw, &rawlen);
    if (rc == 0) return -1;
    if (rc < 0) return 0;

    n = rawlen / 4;
    if (n == 0) return 0;
    /* As for line spectra: length <= 0 is a size query and writes nothing. */
    m = (length > 0 && (size_t)length < n) ? (size_t)length
      : (length > 0 ? n : 0);
    if (profile && m) {
        const unsigned char *p0 = raw, *p1 = raw + n, *p2 = raw + 2 * n,
                            *p3 = raw + 3 * n;
        for (i = 0; i < m; i++)
            profile[i] = (uint32_t)p0[i] | ((uint32_t)p1[i] << 8)
                       | ((uint32_t)p2[i] << 16) | ((uint32_t)p3[i] << 24);
    }
    return (int32_t)n;
}

/* ------------------------------------------------------------- conversion */

uint32_t tsf_index_to_mz(uint64_t handle, int64_t fid, const double *in,
                         double *out, uint32_t cnt) {
    Instance *inst = (Instance *)(uintptr_t)handle;
    const FittedCalib *fc;
    uint32_t i;
    (void)fid;   /* one fit covers the whole analysis */
    g_error[0] = '\0';
    if (!inst) { set_error("Invalid handle."); return 0; }
    if (!in || !out) { set_error("NULL array argument."); return 0; }
    fc = &inst->fitted;
    for (i = 0; i < cnt; i++)
        out[i] = fitted_mz_from_time(fc, fc->delay + in[i] * fc->timebase);
    return 1;
}

uint32_t tsf_mz_to_index(uint64_t handle, int64_t fid, const double *in,
                         double *out, uint32_t cnt) {
    Instance *inst = (Instance *)(uintptr_t)handle;
    const FittedCalib *fc;
    uint32_t i;
    (void)fid;
    g_error[0] = '\0';
    if (!inst) { set_error("Invalid handle."); return 0; }
    if (!in || !out) { set_error("NULL array argument."); return 0; }
    fc = &inst->fitted;
    for (i = 0; i < cnt; i++)
        out[i] = (in[i] >= 0.0)
               ? (fitted_forward(fc->p, in[i]) - fc->delay) / fc->timebase
               : NAN;
    return 1;
}
