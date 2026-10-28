"""Ground-truth oracle: the real Bruker libtimsdata.so, used ONLY for testing.

This is how `pytsf` is validated. On a Linux x86-64 machine the SDK's own shared
library can be loaded and asked the same questions as the pure-Python reader,
and the two answers compared. Nothing in the `pytsf` package itself touches
this file or the Bruker binary -- it is a test fixture.

Set TIMSDATA_SO to the path of libtimsdata.so (or drop it next to this file).
"""

from __future__ import annotations

import ctypes
import os

import numpy as np

_CANDIDATES = [
    os.environ.get("TIMSDATA_SO", ""),
    os.path.join(os.path.dirname(__file__), "libtimsdata.so"),
    "libtimsdata.so",
]


class OracleUnavailable(RuntimeError):
    pass


def _load() -> ctypes.CDLL:
    errors = []
    for cand in _CANDIDATES:
        if not cand:
            continue
        try:
            return ctypes.cdll.LoadLibrary(cand)
        except OSError as exc:
            errors.append(f"{cand}: {exc}")
    raise OracleUnavailable(
        "Could not load the Bruker libtimsdata.so (the SDK ships Linux and Windows builds only; this fixture loads the Linux one). "
        "Set TIMSDATA_SO. Tried:\n  " + "\n  ".join(errors)
    )


class Oracle:
    """Thin ctypes wrapper over the TSF portion of the Bruker C API."""

    def __init__(self) -> None:
        d = _load()
        self.d = d
        d.tsf_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tsf_open.restype = ctypes.c_uint64
        d.tsf_close.argtypes = [ctypes.c_uint64]
        d.tsf_close.restype = None
        d.tsf_get_last_error_string.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tsf_get_last_error_string.restype = ctypes.c_uint32
        d.tsf_has_recalibrated_state.argtypes = [ctypes.c_uint64]
        d.tsf_has_recalibrated_state.restype = ctypes.c_uint32

        conv = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_uint32,
        ]
        d.tsf_index_to_mz.argtypes = conv
        d.tsf_index_to_mz.restype = ctypes.c_uint32
        d.tsf_mz_to_index.argtypes = conv
        d.tsf_mz_to_index.restype = ctypes.c_uint32

        d.tsf_read_line_spectrum_v2.argtypes = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int32,
        ]
        d.tsf_read_line_spectrum_v2.restype = ctypes.c_int32
        d.tsf_read_line_spectrum_with_width_v2.argtypes = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_float),
            ctypes.POINTER(ctypes.c_float),
            ctypes.c_int32,
        ]
        d.tsf_read_line_spectrum_with_width_v2.restype = ctypes.c_int32
        d.tsf_read_profile_spectrum_v2.argtypes = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_int32,
        ]
        d.tsf_read_profile_spectrum_v2.restype = ctypes.c_int32
        self.handle: int | None = None

    # -- lifecycle -------------------------------------------------------

    def last_error(self) -> str:
        buf = ctypes.create_string_buffer(2048)
        self.d.tsf_get_last_error_string(buf, 2048)
        return buf.value.decode(errors="replace")

    def open(self, directory: str, use_recalibrated_state: bool = False) -> int:
        h = self.d.tsf_open(str(directory).encode("utf-8"), 1 if use_recalibrated_state else 0)
        if h == 0:
            raise RuntimeError(f"tsf_open failed: {self.last_error()}")
        self.handle = h
        return h

    def close(self) -> None:
        if self.handle:
            self.d.tsf_close(self.handle)
            self.handle = None

    def __enter__(self) -> Oracle:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # -- queries ---------------------------------------------------------

    def _conv(self, func, frame_id: int, values) -> np.ndarray:
        a = np.ascontiguousarray(values, dtype=np.float64)
        out = np.empty(a.size, dtype=np.float64)
        ok = func(
            self.handle,
            frame_id,
            a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            a.size,
        )
        if not ok:
            raise RuntimeError(self.last_error())
        return out

    def index_to_mz(self, frame_id: int, indices) -> np.ndarray:
        return self._conv(self.d.tsf_index_to_mz, frame_id, indices)

    def mz_to_index(self, frame_id: int, mzs) -> np.ndarray:
        return self._conv(self.d.tsf_mz_to_index, frame_id, mzs)

    def read_line_spectrum(self, frame_id: int, cap: int = 1 << 20):
        idx = np.empty(cap, dtype=np.float64)
        inten = np.empty(cap, dtype=np.float32)
        n = self.d.tsf_read_line_spectrum_v2(
            self.handle,
            frame_id,
            idx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            inten.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            cap,
        )
        if n < 0:
            raise RuntimeError(self.last_error())
        if n > cap:
            return self.read_line_spectrum(frame_id, cap=n)
        return idx[:n].copy(), inten[:n].copy()

    def read_line_spectrum_with_width(self, frame_id: int, cap: int = 1 << 20):
        idx = np.empty(cap, dtype=np.float64)
        inten = np.empty(cap, dtype=np.float32)
        width = np.empty(cap, dtype=np.float32)
        n = self.d.tsf_read_line_spectrum_with_width_v2(
            self.handle,
            frame_id,
            idx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            inten.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            width.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            cap,
        )
        if n < 0:
            raise RuntimeError(self.last_error())
        if n > cap:
            return self.read_line_spectrum_with_width(frame_id, cap=n)
        return idx[:n].copy(), inten[:n].copy(), width[:n].copy()

    def read_profile_spectrum(self, frame_id: int, cap: int = 1 << 21):
        buf = np.empty(cap, dtype=np.uint32)
        n = self.d.tsf_read_profile_spectrum_v2(
            self.handle,
            frame_id,
            buf.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)),
            cap,
        )
        if n < 0:
            raise RuntimeError(self.last_error())
        if n > cap:
            return self.read_profile_spectrum(frame_id, cap=n)
        return buf[:n].copy()


class TdfOracle:
    """Thin ctypes wrapper over the TDF (``tims_*``) portion of the public C API.

    Signatures are those documented in the SDK's ``timsdata.h`` / example
    ``timsdata.py``; this class only *calls* the library.
    """

    def __init__(self) -> None:
        d = _load()
        self.d = d
        d.tims_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tims_open.restype = ctypes.c_uint64
        d.tims_close.argtypes = [ctypes.c_uint64]
        d.tims_close.restype = None
        d.tims_get_last_error_string.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
        d.tims_get_last_error_string.restype = ctypes.c_uint32
        d.tims_read_scans_v2.argtypes = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.c_uint32,
            ctypes.c_uint32,
            ctypes.c_void_p,
            ctypes.c_uint32,
        ]
        d.tims_read_scans_v2.restype = ctypes.c_uint32
        conv = [
            ctypes.c_uint64,
            ctypes.c_int64,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_uint32,
        ]
        for fn in (
            "tims_index_to_mz",
            "tims_mz_to_index",
            "tims_scannum_to_oneoverk0",
            "tims_oneoverk0_to_scannum",
        ):
            getattr(d, fn).argtypes = conv
            getattr(d, fn).restype = ctypes.c_uint32
        d.tims_oneoverk0_to_ccs_for_mz.argtypes = [
            ctypes.c_double,
            ctypes.c_int32,
            ctypes.c_double,
        ]
        d.tims_oneoverk0_to_ccs_for_mz.restype = ctypes.c_double
        self.handle = 0

    def last_error(self) -> str:
        n = self.d.tims_get_last_error_string(None, 0)
        buf = ctypes.create_string_buffer(n)
        self.d.tims_get_last_error_string(buf, n)
        return buf.value.decode(errors="replace")

    def open(self, directory: str, use_recalibrated_state: bool = False) -> int:
        self.handle = self.d.tims_open(directory.encode(), 1 if use_recalibrated_state else 0)
        if self.handle == 0:
            raise RuntimeError(f"tims_open failed: {self.last_error()}")
        return self.handle

    def close(self) -> None:
        if self.handle:
            self.d.tims_close(self.handle)
            self.handle = 0

    def __enter__(self) -> TdfOracle:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def _conv(self, func, frame_id: int, values) -> np.ndarray:
        values = np.ascontiguousarray(values, dtype=np.float64)
        out = np.empty_like(values)
        ok = func(
            self.handle,
            frame_id,
            values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            values.size,
        )
        if not ok:
            raise RuntimeError(self.last_error())
        return out

    def index_to_mz(self, frame_id, indices):
        return self._conv(self.d.tims_index_to_mz, frame_id, indices)

    def mz_to_index(self, frame_id, mzs):
        return self._conv(self.d.tims_mz_to_index, frame_id, mzs)

    def scannum_to_oneoverk0(self, frame_id, scans):
        return self._conv(self.d.tims_scannum_to_oneoverk0, frame_id, scans)

    def oneoverk0_to_scannum(self, frame_id, k0):
        return self._conv(self.d.tims_oneoverk0_to_scannum, frame_id, k0)

    def oneoverk0_to_ccs(self, k0: float, charge: int, mz: float) -> float:
        return self.d.tims_oneoverk0_to_ccs_for_mz(k0, charge, mz)

    def read_scans(self, frame_id: int, scan_begin: int, scan_end: int):
        """Return (scan_of_peak, tof_index, intensity) arrays for a scan range."""
        cnt = 1 << 16
        while True:
            buf = np.empty(cnt, dtype=np.uint32)
            need = self.d.tims_read_scans_v2(self.handle, frame_id, scan_begin, scan_end, buf.ctypes.data, 4 * cnt)
            if need == 0:
                raise RuntimeError(self.last_error())
            if need > 4 * cnt:
                cnt = need // 4 + 1
                continue
            break
        d = scan_end - scan_begin
        counts = buf[:d].astype(np.int64)
        scans, idx, inten = [], [], []
        for i in range(d):
            n = int(counts[i])
            idx.append(buf[d : d + n])
            d += n
            inten.append(buf[d : d + n])
            d += n
            scans.append(np.full(n, scan_begin + i, dtype=np.int32))
        cat = lambda xs, dt: np.concatenate(xs).astype(dt) if xs else np.empty(0, dt)  # noqa: E731
        return cat(scans, np.int32), cat(idx, np.uint32), cat(inten, np.uint32)
