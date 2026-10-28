#!/usr/bin/env python3
"""Differentially validate pytsf's TDF reader against Bruker's library on a
REAL ion-mobility .d dataset.

Needs the vendor library (Linux x86-64 or Windows); it is only *called*.

    TIMSDATA_SO=/path/to/libtimsdata.so \
        python tests/validate_real_tdf.py /path/to/analysis.d [max_frames]

Reports, for every frame: scan structure, TOF indices (must be identical),
intensities (identical, or a constant vendor scale factor -- see
pytsf/tdf.py), m/z axis (ppm, inside the calibrant range), and the mobility
axis (percent, inside the span of the reference ions the file records).
Exits non-zero if the frame decoding disagrees or the m/z axis is off by more
than the tolerance inside the calibrant range.
"""

from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pytsf import TdfData  # noqa: E402

from oracle import TdfOracle  # noqa: E402

MZ_REL_TOL = 5e-6  # 5 ppm inside the calibrant range


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    d = argv[1]
    limit = int(argv[2]) if len(argv) > 2 else 0
    t0 = time.time()
    bad = {"structure": 0, "tof": 0, "intensity": 0}
    n_peaks = n_int_diff = 0
    ratios = []

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tdf = TdfData(d)
    for w in caught:
        print(f"  note: {w.message}")

    with tdf, TdfOracle() as oracle:
        oracle.open(d)
        print(f"{d}")
        print(f"  frames / scans     : {len(tdf)} / {tdf.meta.num_scans}")
        print(
            f"  calibrant range    : {tdf.calibration.mz_min:.1f} - {tdf.calibration.mz_max:.1f} m/z "
            f"({tdf.calibration.ref_mz.size} points, fit RMS {tdf.calibration.fit_rms_ppm:.3f} ppm)"
        )
        print(f"  mobility           : {tdf.mobility}")
        for ion in tdf.reference_ions:
            print(
                f"    {ion.name:>18s}  found={ion.found!s:5s}  S/N={ion.snr:8.1f}  "
                f"m/z err={ion.mz_error_ppm:+6.2f} ppm  apex scan={ion.apex_scan:8.2f}"
            )

        frames = tdf.frames if not limit else tdf.frames[:limit]
        for fr in frames:
            ours = tdf.read_frame(fr.id)
            vs, vi, vint = oracle.read_scans(fr.id, 0, fr.num_scans)
            n_peaks += ours.num_peaks
            if ours.num_peaks != vi.size or not np.array_equal(ours.scan, vs):
                bad["structure"] += 1
                continue
            if not np.array_equal(ours.tof_indices, vi):
                bad["tof"] += 1
            if not np.array_equal(ours.intensities, vint):
                bad["intensity"] += 1
                n_int_diff += int((ours.intensities != vint).sum())
                m = ours.intensities > 0
                if m.any():
                    ratios.append(np.median(vint[m] / ours.intensities[m]))

        # m/z axis
        fid = frames[0].id
        idx = np.linspace(0, (tdf.meta.digitizer_num_samples or 400000) - 1, 4001)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ours_mz = tdf.index_to_mz(idx)
        theirs_mz = oracle.index_to_mz(fid, idx)
        lo, hi = tdf.calibration.mz_min, tdf.calibration.mz_max
        inside = (theirs_mz >= lo) & (theirs_mz <= hi) & np.isfinite(ours_mz)
        ppm = np.abs(ours_mz - theirs_mz) / theirs_mz * 1e6
        if inside.any():
            print(f"  m/z inside range   : median {np.median(ppm[inside]):.3f} ppm  max {ppm[inside].max():.3f} ppm")
        else:
            print(
                "  m/z                : the vendor library places no index inside the calibrant range "
                "(a synthetic file with placeholder MzCalibration?); m/z not compared"
            )
        if inside.any() and (~inside).any():
            print(
                f"  m/z extrapolated   : median {np.nanmedian(ppm[~inside]):.1f} ppm  "
                f"max {np.nanmax(ppm[~inside]):.1f} ppm  (outside the calibrant range; expected)"
            )

        # mobility axis
        scans = np.arange(tdf.meta.num_scans, dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ours_k0 = tdf.scan_to_oneoverk0(scans)
        try:
            theirs_k0 = oracle.scannum_to_oneoverk0(fid, scans)
        except RuntimeError as exc:
            print(f"  1/K0               : vendor library cannot convert for this file ({exc}); not compared")
        else:
            pct = np.abs(ours_k0 - theirs_k0) / theirs_k0 * 100
            mob = tdf.mobility
            span = (scans >= mob.scan_min) & (scans <= mob.scan_max)
            tag = "reference ions" if mob.method == "reference_ions" else "acq-range approximation"
            print(
                f"  1/K0 ({tag:>22s}): median {np.median(pct[span]):.3f}%  max {pct[span].max():.3f}%  "
                f"(scans {int(mob.scan_min)}-{int(mob.scan_max)})"
            )

    print(f"  peaks compared     : {n_peaks:,} in {len(frames)} frames, {time.time() - t0:.1f} s")
    ok = True
    if bad["structure"] or bad["tof"]:
        print(f"  FAIL: scan structure/TOF mismatch in {bad['structure'] + bad['tof']} frames")
        ok = False
    if bad["intensity"]:
        r = np.median(ratios) if ratios else float("nan")
        same = ratios and np.allclose(ratios, r, rtol=1e-3)
        print(
            f"  intensities        : vendor = stored x {r:.4g}; {n_int_diff:,} peaks differ in "
            f"{bad['intensity']} frames ({'constant factor' if same else 'NOT a constant factor'}); "
            f"pytsf returns the stored values"
        )
        if not same:
            ok = False
    else:
        print("  intensities        : identical")
    if inside.sum() and ppm[inside].max() > MZ_REL_TOL * 1e6:
        print("  FAIL: m/z disagrees inside the calibrant range")
        ok = False
    print(
        "  RESULT:",
        "PASS -- frames decode identically; m/z agrees within tolerance" if ok else "FAIL",
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
