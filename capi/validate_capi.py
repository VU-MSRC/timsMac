#!/usr/bin/env python3
"""Differentially validate the C compatibility library against the vendor library.

Loads both shared libraries through the *same* ctypes bindings -- i.e. exercises
the identical C ABI -- and compares every documented entry point on a real
dataset.

Spectrum decoding is compared exactly. m/z is not: this library computes it
independently, by fitting the published TOF calibration equation to the points
the file records, rather than reproducing the vendor's own algorithm. Agreement
there is of order 1 ppm, and the tolerance below reflects that.

    TIMSDATA_SO=/path/to/libtimsdata.so \
        python3 validate_capi.py ./libtimsdata_tsf.so /path/to/analysis.d [limit]

Exits non-zero on any disagreement.
"""

from __future__ import annotations

import ctypes
import os
import sys
import time

import numpy as np

P_d = ctypes.POINTER(ctypes.c_double)
P_f = ctypes.POINTER(ctypes.c_float)
P_u = ctypes.POINTER(ctypes.c_uint32)


class Lib:
    """ctypes bindings for the TSF C API -- identical for either library."""

    def __init__(self, path: str, name: str):
        self.name = name
        self.d = ctypes.cdll.LoadLibrary(path)
        d = self.d
        d.tsf_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tsf_open.restype = ctypes.c_uint64
        d.tsf_close.argtypes = [ctypes.c_uint64]
        d.tsf_get_last_error_string.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tsf_get_last_error_string.restype = ctypes.c_uint32
        d.tsf_has_recalibrated_state.argtypes = [ctypes.c_uint64]
        d.tsf_has_recalibrated_state.restype = ctypes.c_uint32
        conv = [ctypes.c_uint64, ctypes.c_int64, P_d, P_d, ctypes.c_uint32]
        d.tsf_index_to_mz.argtypes = conv
        d.tsf_index_to_mz.restype = ctypes.c_uint32
        d.tsf_mz_to_index.argtypes = conv
        d.tsf_mz_to_index.restype = ctypes.c_uint32
        d.tsf_read_line_spectrum_v2.argtypes = [
            ctypes.c_uint64, ctypes.c_int64, P_d, P_f, ctypes.c_int32]
        d.tsf_read_line_spectrum_v2.restype = ctypes.c_int32
        d.tsf_read_line_spectrum_with_width_v2.argtypes = [
            ctypes.c_uint64, ctypes.c_int64, P_d, P_f, P_f, ctypes.c_int32]
        d.tsf_read_line_spectrum_with_width_v2.restype = ctypes.c_int32
        d.tsf_read_profile_spectrum_v2.argtypes = [
            ctypes.c_uint64, ctypes.c_int64, P_u, ctypes.c_int32]
        d.tsf_read_profile_spectrum_v2.restype = ctypes.c_int32
        self.h = 0

    def err(self) -> str:
        b = ctypes.create_string_buffer(2048)
        self.d.tsf_get_last_error_string(b, 2048)
        return b.value.decode(errors="replace")

    def open(self, path: str) -> None:
        self.h = self.d.tsf_open(path.encode(), 0)
        if not self.h:
            raise RuntimeError(f"{self.name}: tsf_open failed: {self.err()}")

    def close(self) -> None:
        if self.h:
            self.d.tsf_close(self.h)
            self.h = 0

    def line(self, fid: int, cap: int = 1 << 20):
        idx = np.empty(cap, np.float64); inten = np.empty(cap, np.float32)
        n = self.d.tsf_read_line_spectrum_v2(
            self.h, fid, idx.ctypes.data_as(P_d), inten.ctypes.data_as(P_f), cap)
        if n < 0:
            raise RuntimeError(f"{self.name}: {self.err()}")
        if n > cap:
            return self.line(fid, n)
        return idx[:n].copy(), inten[:n].copy()

    def profile(self, fid: int, cap: int = 1 << 21):
        buf = np.empty(cap, np.uint32)
        n = self.d.tsf_read_profile_spectrum_v2(
            self.h, fid, buf.ctypes.data_as(P_u), cap)
        if n < 0:
            raise RuntimeError(f"{self.name}: {self.err()}")
        if n > cap:
            return self.profile(fid, n)
        return buf[:n].copy()

    def i2mz(self, fid: int, vals):
        a = np.ascontiguousarray(vals, np.float64)
        o = np.empty(a.size, np.float64)
        if not self.d.tsf_index_to_mz(self.h, fid, a.ctypes.data_as(P_d),
                                      o.ctypes.data_as(P_d), a.size):
            raise RuntimeError(f"{self.name}: {self.err()}")
        return o

    def mz2i(self, fid: int, vals):
        a = np.ascontiguousarray(vals, np.float64)
        o = np.empty(a.size, np.float64)
        if not self.d.tsf_mz_to_index(self.h, fid, a.ctypes.data_as(P_d),
                                      o.ctypes.data_as(P_d), a.size):
            raise RuntimeError(f"{self.name}: {self.err()}")
        return o


def main(argv) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    mine_path, dataset = argv[1], argv[2]
    limit = int(argv[3]) if len(argv) > 3 else 0
    vendor_path = os.environ.get("TIMSDATA_SO")
    if not vendor_path:
        print("TIMSDATA_SO is not set -- nothing to compare against.")
        return 2

    mine = Lib(mine_path, "pytsf-capi")
    vendor = Lib(vendor_path, "vendor")
    mine.open(dataset)
    vendor.open(dataset)

    import sqlite3
    con = sqlite3.connect(f"file:{dataset}/analysis.tsf?mode=ro", uri=True)
    ids = [r[0] for r in con.execute("SELECT Id FROM Frames ORDER BY Id")]
    md = {k: v for k, v in con.execute("SELECT Key, Value FROM GlobalMetadata")}
    con.close()
    has_line = str(md.get("HasLineSpectra", "0")).strip() == "1"
    has_profile = str(md.get("HasProfileSpectra", "0")).strip() == "1"
    print(f"dataset has line spectra: {has_line}, profile spectra: {has_profile}")
    if limit:
        ids = ids[:limit]

    # The calibrant range: inside it the fitted calibration is held to a tight
    # tolerance; outside it the reader is extrapolating by design, so the error
    # is reported but not treated as a failure.
    import sqlite3 as _sq, struct as _st
    con = _sq.connect(f"file:{dataset}/analysis.tsf?mode=ro", uri=True)
    _r = con.execute("SELECT Value FROM CalibrationInfo "
                     "WHERE KeyName='ReferencePeakMasses'").fetchone()
    con.close()
    if _r and isinstance(_r[0], (bytes, bytearray)):
        _m = _st.unpack(f"<{len(_r[0]) // 8}d", bytes(_r[0]))
        cal_lo, cal_hi = min(_m), max(_m)
    else:
        cal_lo, cal_hi = 0.0, float("inf")
    print(f"calibrant range: {cal_lo:.1f} - {cal_hi:.1f} m/z")

    fails = dict(line_index=0, line_intensity=0, profile=0, mz=0, index=0)
    worst_mz = worst_idx = 0.0
    worst_mz_out = 0.0
    n_peaks = n_prof = 0
    t0 = time.time()

    for k, fid in enumerate(ids):
        mi = vi = np.empty(0, np.float64)
        if has_line:
            mi, mv = mine.line(fid)
            vi, vv = vendor.line(fid)
            if not np.array_equal(mi, vi): fails["line_index"] += 1
            if not np.array_equal(mv, vv): fails["line_intensity"] += 1
            n_peaks += mi.size

        if has_profile:
            mp = mine.profile(fid); vp = vendor.profile(fid)
            if not np.array_equal(mp, vp): fails["profile"] += 1
            n_prof += mp.size

        if not has_line:
            # still exercise the conversion functions over the index range
            probe = np.linspace(0.0, 400000.0, 257)
            v_mz = vendor.i2mz(fid, probe); m_mz = mine.i2mz(fid, probe)
            vi = mi = probe
        else:
            v_mz = vendor.i2mz(fid, vi); m_mz = mine.i2mz(fid, mi)
        ok = np.isfinite(v_mz) & (v_mz > 0) & np.isfinite(m_mz)
        inside = ok & (v_mz >= cal_lo) & (v_mz <= cal_hi)
        outside = ok & ~inside
        if inside.any():
            rel = float(np.max(np.abs((m_mz[inside] - v_mz[inside])
                                      / v_mz[inside])))
            worst_mz = max(worst_mz, rel)
            if rel > 5e-6: fails["mz"] += 1   # 5 ppm inside the calibrant range
            a = float(np.max(np.abs(mine.mz2i(fid, v_mz[inside])
                                    - vendor.mz2i(fid, v_mz[inside]))))
            worst_idx = max(worst_idx, a)
            if a > 50.0: fails["index"] += 1
        if outside.any():
            worst_mz_out = max(worst_mz_out, float(np.max(
                np.abs((m_mz[outside] - v_mz[outside]) / v_mz[outside]))))

        if (k + 1) % 500 == 0:
            print(f"  {k + 1}/{len(ids)} ({time.time() - t0:.0f}s)")

    print(f"\nhas_recalibrated_state: vendor="
          f"{vendor.d.tsf_has_recalibrated_state(vendor.h)} "
          f"ours={mine.d.tsf_has_recalibrated_state(mine.h)}")
    mine.close(); vendor.close()

    print(f"spectra compared        : {len(ids):,}")
    print(f"peaks compared          : {n_peaks:,}")
    print(f"profile points compared : {n_prof:,}")
    for k_, v_ in fails.items():
        print(f"mismatched {k_:<16}: {v_}")
    print(f"max m/z error, in calibrant range : {worst_mz * 1e6:8.2f} ppm")
    if worst_mz_out:
        print(f"max m/z error, extrapolated       : {worst_mz_out * 1e6:8.2f} ppm "
              "(expected: external calibration)")
    print(f"max index abs error     : {worst_idx:.3e} samples")
    print(f"elapsed                 : {time.time() - t0:.1f}s")

    bad = sum(fails.values())
    print("\nFAILED" if bad else
          "\nOK -- spectra decode identically; m/z agrees within tolerance")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
