"""Exercise the installed C API without requiring the vendor library."""

# ruff: file-ignore[assert]

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

import numpy as np
import pytest
import pytsf
from pytsf import TsfData

import synth


@pytest.fixture(scope="module")
def capi() -> ctypes.CDLL:
    """Load the packaged shared library, requiring it in native wheel tests."""
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    library = Path(pytsf.__file__).with_name(f"libtimsdata_tsf{suffix}")
    if not library.exists():
        if os.environ.get("PYTSF_REQUIRE_CAPI") == "1":
            pytest.fail(f"The installed wheel contains no C API library: {library}")
        pytest.skip("The C API is only included in native wheels.")
    assert library.with_name("tsf_compat.h").is_file()
    lib = ctypes.CDLL(str(library))
    lib.tsf_compat_version.argtypes = []
    lib.tsf_compat_version.restype = ctypes.c_char_p
    lib.tsf_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
    lib.tsf_open.restype = ctypes.c_uint64
    lib.tsf_close.argtypes = [ctypes.c_uint64]
    lib.tsf_close.restype = None
    lib.tsf_read_line_spectrum_v2.argtypes = [
        ctypes.c_uint64,
        ctypes.c_int64,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_float),
        ctypes.c_int32,
    ]
    lib.tsf_read_line_spectrum_v2.restype = ctypes.c_int32
    lib.tsf_index_to_mz.argtypes = [
        ctypes.c_uint64,
        ctypes.c_int64,
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_uint32,
    ]
    lib.tsf_index_to_mz.restype = ctypes.c_uint32
    return lib


def test_capi_version(capi: ctypes.CDLL) -> None:
    """The wheel exposes the standalone C ABI and its version entry point."""
    assert capi.tsf_compat_version().startswith(b"pytsf-capi ")


def test_capi_spectrum(capi: ctypes.CDLL, tmp_path: Path) -> None:
    """Compare C spectrum decoding and calibration with the Python reader."""
    calibration = synth.Calibration()
    indices = calibration.mz_to_index(np.array([500.0, 900.0, 1200.0]))
    dataset = synth.build(
        str(tmp_path / "spectrum.d"),
        [synth.Spectrum(indices.tolist(), [100.0, 250.0, 1000.0])],
        calibration=calibration,
        has_profile=False,
    )
    handle = capi.tsf_open(dataset.encode(), 0)
    assert handle != 0
    try:
        count = capi.tsf_read_line_spectrum_v2(handle, 1, None, None, 0)
        assert count == 3
        index = (ctypes.c_double * count)()
        intensity = (ctypes.c_float * count)()
        mz = (ctypes.c_double * count)()
        assert capi.tsf_read_line_spectrum_v2(handle, 1, index, intensity, count) == count
        assert capi.tsf_index_to_mz(handle, 1, index, mz, count) == 1
        with TsfData(dataset) as tsf:
            expected_index, expected_intensity = tsf.read_line_spectrum(1)
            np.testing.assert_array_equal(np.ctypeslib.as_array(index), expected_index)
            np.testing.assert_array_equal(np.ctypeslib.as_array(intensity), expected_intensity)
            np.testing.assert_allclose(np.ctypeslib.as_array(mz), tsf.index_to_mz(1, expected_index), rtol=1e-6)
    finally:
        capi.tsf_close(handle)


def test_capi_rejects_missing_dataset(capi: ctypes.CDLL, tmp_path: Path) -> None:
    """A failed open reports failure through the C ABI."""
    assert capi.tsf_open(str(tmp_path / "missing.d").encode(), 0) == 0
