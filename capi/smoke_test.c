/* smoke_test.c -- minimal check that the compatibility library works.
 *
 *   ./smoke_test /path/to/analysis.d
 *
 * Opens the dataset, reads the first spectrum, converts indices to m/z, and
 * prints a short summary. Exit status 0 means the library is working.
 *
 * Also serves as a worked example of calling the library from C.
 */

#include <stdio.h>
#include <stdlib.h>

#include "tsf_compat.h"

static void die(const char *what) {
    char buf[1024];
    tsf_get_last_error_string(buf, sizeof(buf));
    fprintf(stderr, "%s: %s\n", what, buf);
    exit(1);
}

int main(int argc, char **argv) {
    uint64_t h;
    int32_t n, got;
    double *idx = NULL, *mz = NULL;
    float *inten = NULL;
    uint32_t *profile = NULL;
    int i;

    if (argc < 2) {
        fprintf(stderr, "usage: %s /path/to/analysis.d\n", argv[0]);
        return 2;
    }

    printf("%s\n", tsf_compat_version());

    h = tsf_open(argv[1], 0);
    if (!h) die("tsf_open");
    printf("opened %s\n", argv[1]);

    /* --- line spectrum ------------------------------------------------- */
    n = tsf_read_line_spectrum_v2(h, 1, NULL, NULL, 0);
    if (n < 0) die("tsf_read_line_spectrum_v2 (sizing)");
    printf("spectrum 1: %d peaks\n", n);

    if (n > 0) {
        idx   = (double *)malloc((size_t)n * sizeof(double));
        mz    = (double *)malloc((size_t)n * sizeof(double));
        inten = (float  *)malloc((size_t)n * sizeof(float));
        if (!idx || !mz || !inten) { fprintf(stderr, "out of memory\n"); return 1; }

        got = tsf_read_line_spectrum_v2(h, 1, idx, inten, n);
        if (got < 0) die("tsf_read_line_spectrum_v2");

        if (!tsf_index_to_mz(h, 1, idx, mz, (uint32_t)got))
            die("tsf_index_to_mz");

        printf("  m/z range : %.4f .. %.4f\n", mz[0], mz[got - 1]);
        printf("  first 5 peaks:\n");
        for (i = 0; i < (got < 5 ? got : 5); i++)
            printf("    mz=%14.6f  intensity=%12.1f\n", mz[i], (double)inten[i]);

        /* round-trip the m/z back to indices */
        {
            double *back = (double *)malloc((size_t)got * sizeof(double));
            double worst = 0.0;
            if (back) {
                if (!tsf_mz_to_index(h, 1, mz, back, (uint32_t)got))
                    die("tsf_mz_to_index");
                for (i = 0; i < got; i++) {
                    double d = back[i] - idx[i];
                    if (d < 0) d = -d;
                    if (d > worst) worst = d;
                }
                printf("  index -> m/z -> index, worst drift: %.3e samples\n", worst);
                free(back);
            }
        }
    }

    /* --- profile spectrum ---------------------------------------------- */
    n = tsf_read_profile_spectrum_v2(h, 1, NULL, 0);
    if (n < 0) {
        char buf[1024];
        tsf_get_last_error_string(buf, sizeof(buf));
        printf("profile: not available (%s)\n", buf);
    } else {
        printf("profile: %d points\n", n);
        if (n > 0) {
            long long sum = 0;
            uint32_t peak = 0;
            profile = (uint32_t *)malloc((size_t)n * sizeof(uint32_t));
            if (profile) {
                if (tsf_read_profile_spectrum_v2(h, 1, profile, n) < 0)
                    die("tsf_read_profile_spectrum_v2");
                for (i = 0; i < n; i++) {
                    sum += profile[i];
                    if (profile[i] > peak) peak = profile[i];
                }
                printf("  summed intensity %lld, base peak %u\n", sum, peak);
            }
        }
    }

    free(idx); free(mz); free(inten); free(profile);
    tsf_close(h);
    printf("\nOK\n");
    return 0;
}
