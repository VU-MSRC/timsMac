"""Differential test: pytsf vs Bruker's libtimsdata.so on identical input.

Skips automatically when the vendor library is unavailable (e.g. on macOS,
which is the whole reason this package exists). Set TIMSDATA_SO to run it.
"""

from __future__ import annotations

import os
import random
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import synth  # noqa: E402  (same directory)
from pytsf import TsfData  # noqa: E402

try:
    from oracle import Oracle, OracleUnavailable
    try:
        Oracle()
        HAVE_ORACLE = True
        SKIP_REASON = ""
    except OracleUnavailable as exc:
        HAVE_ORACLE = False
        SKIP_REASON = str(exc)
except ImportError as exc:  # pragma: no cover
    HAVE_ORACLE = False
    SKIP_REASON = f"could not import oracle: {exc}"

needs_oracle = pytest.mark.skipif(not HAVE_ORACLE, reason=SKIP_REASON)

INDICES = np.concatenate([
    np.array([0.0, 0.5, 1.0, 123.456, 399999.0, 400000.0]),
    np.linspace(0.0, 400000.0, 257),
])


def _spectra(seed: int = 4242):
    rnd = random.Random(seed)
    out = []
    for n, plen in [(1, 1), (2, 5), (7, 100), (333, 4096), (1024, 20001), (0, 3)]:
        idx = sorted(rnd.uniform(0, 399999) for _ in range(n))
        # integers below 2**24 are exact in float32, so round-trips are exact
        inten = [float(rnd.randint(1, 100000)) for _ in range(n)]
        width = [rnd.uniform(1e-4, 5.0) for _ in range(n)]
        prof = [0] * plen
        for _ in range(min(plen, 40)):
            prof[rnd.randrange(plen)] = rnd.choice(
                [1, 255, 256, 65535, 65536, 16777215, 2 ** 31, 2 ** 32 - 1]
            )
        out.append(synth.Spectrum(idx, inten, width, prof,
                                  gap_before=rnd.randrange(0, 17)))
    return out


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("tsf") / "analysis.d"
    return synth.build(str(d), _spectra())


# -- blob decoding --------------------------------------------------------

@needs_oracle
def test_line_spectra_match_vendor(dataset):
    with TsfData(dataset) as tsf, Oracle() as oracle:
        oracle.open(dataset)
        for frame in tsf.frames:
            ours = tsf.read_line_spectrum_with_width(frame.id)
            theirs = oracle.read_line_spectrum_with_width(frame.id)
            for name, a, b in zip(("index", "intensity", "width"), ours, theirs):
                assert np.array_equal(a, b), f"{name} differs on spectrum {frame.id}"


@needs_oracle
def test_profile_spectra_match_vendor(dataset):
    with TsfData(dataset) as tsf, Oracle() as oracle:
        oracle.open(dataset)
        for frame in tsf.frames:
            ours = tsf.read_profile_spectrum(frame.id)
            theirs = oracle.read_profile_spectrum(frame.id)
            assert np.array_equal(ours, theirs), f"profile differs on {frame.id}"


# -- calibration ----------------------------------------------------------

CAL_FAMILIES = [
    ("linear",    synth.Calibration(ml1=1.5e6, ml2=1234.5, ml3=0.0)),
    ("quadratic", synth.Calibration(ml1=1.5e6, ml2=1234.5, ml3=-2e-4)),
    ("shifted",   synth.Calibration(ml1=9.0e5, ml2=-500.0, ml3=1e-4)),
]


@pytest.mark.parametrize("name,cal", CAL_FAMILIES, ids=[n for n, _ in CAL_FAMILIES])
def test_fitted_calibration_recovers_the_published_equation(tmp_path, name, cal):
    """Data synthesised from the published equation must be recovered by the
    fit to numerical precision -- the fit and the generator are the same model."""
    spectra = [synth.Spectrum([1.0], [1.0], [0.1], [0])]
    d = synth.build(str(tmp_path / f"{name}.d"), spectra, calibration=cal)
    with TsfData(d) as tsf:
        probe_mz = np.linspace(400.0, 1700.0, 300)
        idx = cal.mz_to_index(probe_mz)
        got = tsf.index_to_mz(1, idx)
        ppm = np.abs((got - probe_mz) / probe_mz) * 1e6
        assert ppm.max() < 1e-3, f"{name}: max {ppm.max():.5f} ppm"


@pytest.mark.filterwarnings("ignore:.*outside the calibrant range:RuntimeWarning")
@pytest.mark.parametrize("name,cal", CAL_FAMILIES, ids=[n for n, _ in CAL_FAMILIES])
def test_fitted_calibration_round_trips(tmp_path, name, cal):
    spectra = [synth.Spectrum([1.0], [1.0], [0.1], [0])]
    d = synth.build(str(tmp_path / f"{name}rt.d"), spectra, calibration=cal)
    with TsfData(d) as tsf:
        mz = tsf.index_to_mz(1, INDICES)
        sel = np.isfinite(mz) & (mz > 0)
        assert sel.any()
        back = tsf.mz_to_index(1, mz[sel])
        assert np.abs(back - INDICES[sel]).max() < 1e-6


# -- self-consistency (runs without the vendor library) -------------------

@pytest.mark.filterwarnings("ignore:.*outside the calibrant range:RuntimeWarning")
def test_round_trip_without_oracle(dataset):
    """index -> m/z -> index is stable to well under one digitizer sample."""
    with TsfData(dataset) as tsf:
        mz = tsf.index_to_mz(1, INDICES)
        phys = np.isfinite(mz) & (mz > 0)
        assert phys.any()
        back = tsf.mz_to_index(1, mz[phys])
        assert np.abs(back - INDICES[phys]).max() < 1e-6


def test_encode_decode_round_trip(dataset):
    """The decoder recovers exactly what the encoder wrote."""
    spectra = _spectra()
    with TsfData(dataset) as tsf:
        for frame, spec in zip(tsf.frames, spectra):
            idx, inten, width = tsf.read_line_spectrum_with_width(frame.id)
            assert np.array_equal(idx, np.asarray(spec.index, dtype=np.float64))
            assert np.array_equal(inten, np.asarray(spec.intensity, dtype=np.float32))
            assert np.array_equal(width, np.asarray(spec.width, dtype=np.float32))
            prof = tsf.read_profile_spectrum(frame.id)
            assert np.array_equal(prof, np.asarray(spec.profile, dtype=np.uint32))


def test_metadata(dataset):
    with TsfData(dataset) as tsf:
        assert len(tsf) == 6
        assert tsf.meta.has_line_spectra
        assert tsf.meta.has_profile_spectra
        assert tsf.meta.digitizer_num_samples == 400000
        assert tsf.meta.mz_calibration(1).digitizer_timebase > 0


def test_maldi_table_accessible(tmp_path):
    """Files from MALDI imaging runs carry a MaldiFrameInfo table that is not
    part of the documented core schema; it must still be readable."""
    import sqlite3
    spectra = [synth.Spectrum([1.0], [1.0], [0.1], [0])]
    d = synth.build(str(tmp_path / "maldi.d"), spectra)
    con = sqlite3.connect(d + "/analysis.tsf")
    con.execute("CREATE TABLE MaldiFrameInfo (Frame INTEGER, SpotName TEXT, "
                "XIndexPos INTEGER, YIndexPos INTEGER)")
    con.execute("INSERT INTO MaldiFrameInfo VALUES (1,'A1',3,7)")
    con.commit(); con.close()
    with TsfData(d) as tsf:
        assert tsf.meta.is_maldi
        rows = tsf.meta.table("MaldiFrameInfo")
        assert rows[0]["SpotName"] == "A1"
        assert (rows[0]["XIndexPos"], rows[0]["YIndexPos"]) == (3, 7)


# -- the fitted calibration -----------------------------------------------

def test_calibration_is_fitted_from_the_file(dataset):
    with TsfData(dataset) as tsf:
        assert tsf.calibration.ref_mz.size >= 3


def test_fails_loudly_without_calibration_points(tmp_path):
    """With no calibration points in the file there is no way to compute m/z;
    the reader must say so at open time rather than guess."""
    from pytsf.calibration_fitted import CalibrationPointsUnavailable
    d = synth.build(str(tmp_path / "nopoints.d"),
                    [synth.Spectrum([1.0], [1.0], [0.1], [0])],
                    reference_points=False)
    with pytest.raises(CalibrationPointsUnavailable):
        TsfData(d)


def test_fitted_calibration_warns_outside_calibrant_range(tmp_path):
    """Accuracy degrades when extrapolating beyond the calibrant; warn."""
    from pytsf.calibration_fitted import FittedCalibration
    cal = synth.Calibration()
    d = synth.build(str(tmp_path / "warn.d"),
                    [synth.Spectrum([1.0], [1.0], [0.1], [0])], calibration=cal)
    fitted = FittedCalibration(d + "/analysis.tsf")
    with pytest.warns(RuntimeWarning, match="outside the calibrant range"):
        fitted.index_to_mz(cal.mz_to_index(np.array([120.0, 900.0])))
