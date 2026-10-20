#!/usr/bin/env python3
"""Differentially validate pytsf against Bruker's library on a REAL .d dataset.

This is the test that matters before trusting pytsf on data it has not seen.
It needs the vendor library, so run it on Linux x86-64 (or Windows) with the
SDK available -- not on the Mac where you will actually use pytsf.

    TIMSDATA_SO=/path/to/libtimsdata.so \
        python tests/validate_real_dataset.py /path/to/analysis.d

Compares, for every spectrum: line indices, line intensities, peak widths
(when present), profile arrays, index->m/z and m/z->index. Exits non-zero if
anything disagrees.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from oracle import Oracle  # noqa: E402
from pytsf import TsfData  # noqa: E402

# m/z is computed independently from the published equation, not copied from
# the vendor, so agreement is of order 1 ppm rather than exact.
MZ_REL_TOL = 5e-6      # 5 ppm
IDX_ABS_TOL = 50.0     # digitizer samples


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    d = argv[1]
    limit = int(argv[2]) if len(argv) > 2 else 0

    fails = {k: 0 for k in
             ("line_index", "line_intensity", "width", "profile", "mz", "index")}
    worst_mz = 0.0
    worst_idx = 0.0
    worst_out = 0.0
    n_peaks = n_prof = 0
    t0 = time.time()

    with TsfData(d) as tsf, Oracle() as oracle:
        oracle.open(d)
        meta = tsf.meta

        print(f"{d}")
        print(f"  spectra            : {len(tsf)}")
        print(f"  line / profile     : {meta.has_line_spectra} / "
              f"{meta.has_profile_spectra}  (widths: "
              f"{meta.has_line_spectra_peak_width})")
        print(f"  MALDI imaging      : {meta.is_maldi}")
        print(f"  instrument         : "
              f"{meta.global_metadata.get('InstrumentName')!r} "
              f"{meta.global_metadata.get('AcquisitionSoftwareVersion')!r}")
        print(f"  calibration        : {tsf.calibration}")
        print()

        ids = meta.frame_ids()
        if limit:
            ids = ids[:limit]

        for n, fid in enumerate(ids):
            if meta.has_line_spectra:
                if meta.has_line_spectra_peak_width:
                    oi, ov, ow = tsf.read_line_spectrum_with_width(fid)
                    ti, tv, tw = oracle.read_line_spectrum_with_width(fid)
                    if not np.array_equal(ow, tw):
                        fails["width"] += 1
                else:
                    oi, ov = tsf.read_line_spectrum(fid)
                    ti, tv = oracle.read_line_spectrum(fid)
                if not np.array_equal(oi, ti):
                    fails["line_index"] += 1
                if not np.array_equal(ov, tv):
                    fails["line_intensity"] += 1
                n_peaks += oi.size

                tmz = oracle.index_to_mz(fid, ti)
                omz = tsf.index_to_mz(fid, oi)
                ok = np.isfinite(tmz) & (tmz > 0) & np.isfinite(omz)
                lo, hi = tsf.calibration.mz_min, tsf.calibration.mz_max
                m = ok & (tmz >= lo) & (tmz <= hi)
                out = ok & ~m
                if out.any():
                    worst_out = max(worst_out, float(np.max(
                        np.abs((omz[out] - tmz[out]) / tmz[out]))))
                if m.any():
                    rel = float(np.max(np.abs((omz[m] - tmz[m]) / tmz[m])))
                    worst_mz = max(worst_mz, rel)
                    if rel > MZ_REL_TOL:
                        fails["mz"] += 1
                    oidx = tsf.mz_to_index(fid, tmz[m])
                    tidx = oracle.mz_to_index(fid, tmz[m])
                    a = float(np.max(np.abs(oidx - tidx)))
                    worst_idx = max(worst_idx, a)
                    if a > IDX_ABS_TOL:
                        fails["index"] += 1

            if meta.has_profile_spectra:
                op = tsf.read_profile_spectrum(fid)
                tp = oracle.read_profile_spectrum(fid)
                if not np.array_equal(op, tp):
                    fails["profile"] += 1
                n_prof += op.size

            if (n + 1) % 500 == 0:
                print(f"  {n + 1}/{len(ids)} spectra ({time.time() - t0:.0f}s)")

    print()
    print(f"peaks compared          : {n_peaks:,}")
    print(f"profile points compared : {n_prof:,}")
    for k, v in fails.items():
        print(f"mismatched {k:<21}: {v}")
    print(f"max m/z error, in calibrant range : {worst_mz * 1e6:8.2f} ppm")
    if worst_out:
        print(f"max m/z error, extrapolated       : {worst_out * 1e6:8.2f} ppm "
              "(expected: external calibration)")
    print(f"max index abs error     : {worst_idx:.3e} samples")
    print(f"elapsed                 : {time.time() - t0:.1f}s")

    bad = sum(fails.values())
    print()
    print("FAILED" if bad else
          "OK -- spectra decode identically; m/z agrees within tolerance\n"
          "     (computed independently from the published equation, not the\n"
          "      vendor's algorithm -- see PROVENANCE.md)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
