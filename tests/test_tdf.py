"""TDF (ion-mobility) tests: synthetic analyses written by tests/synth.py,
decoded by pytsf and -- when the vendor library is available -- by the vendor
library, and compared.

The vendor library is only ever *called* here; see tests/oracle.py.
"""

from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import synth  # noqa: E402
from pytsf import TdfData, TsfError, MobilityCalibration, oneoverk0_to_ccs, ccs_to_oneoverk0  # noqa: E402
from pytsf.tdf_bin import decode_frame_payload, encode_frame_payload, TdfBinError  # noqa: E402
from test_against_vendor import HAVE_ORACLE, SKIP_REASON  # noqa: E402

if HAVE_ORACLE:
    from oracle import TdfOracle

needs_oracle = pytest.mark.skipif(not HAVE_ORACLE, reason=SKIP_REASON)

NUM_SCANS = 500
# a synthetic "tune mix": (name, m/z, 1/K0) -- values are arbitrary but
# monotone the way a real calibrant is
REFERENCE = [("a", 322.0481, 0.7363), ("b", 622.0290, 0.9915),
             ("c", 922.0098, 1.1986), ("d", 1221.9906, 1.3934), ("e", 1521.9715, 1.5685)]


def _frames(seed=7):
    rnd = np.random.default_rng(seed)
    out = []
    # frame 1: sparse; frame 2: empty scans everywhere; frame 3: busy; frame 4: empty frame
    for n in (50, 5, 4000, 0):
        scan = rnd.integers(0, NUM_SCANS, n)
        mz = rnd.uniform(150, 1650, n)
        inten = rnd.integers(1, 2 ** 20, n)
        out.append(synth.TdfFrameSpec(scan, mz, inten))
    # a PASEF MS/MS frame with a PasefFrameMsMsInfo row
    scan = rnd.integers(100, 300, 30)
    out.append(synth.TdfFrameSpec(scan, rnd.uniform(150, 900, 30), rnd.integers(1, 1000, 30),
                                  msms_type=8, pasef=[(100, 300, 622.03, 2.0, 30.0, 1)]))
    return out


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("tdf") / "analysis.d"
    return synth.build_tdf(str(d), _frames(), num_scans=NUM_SCANS, reference_ions=REFERENCE)


@pytest.fixture(scope="module")
def blank(tmp_path_factory):
    """Same calibrant recorded, but no calibrant ions present in the data."""
    d = tmp_path_factory.mktemp("tdf_blank") / "analysis.d"
    frames = _frames(seed=11)[:3]
    synth.build_tdf(str(d), frames, num_scans=NUM_SCANS, reference_ions=None)
    # record the calibrant in CalibrationInfo without writing its peaks
    import sqlite3, struct
    con = sqlite3.connect(os.path.join(str(d), "analysis.tdf"))
    con.execute("INSERT INTO CalibrationInfo VALUES ('+','ReferenceMobilityPeakNames',?)",
                (b"".join(n.encode() + b"\x00" for n, _, _ in REFERENCE),))
    con.execute("INSERT INTO CalibrationInfo VALUES ('+','ReferencePeakMobilities',?)",
                (struct.pack("<5d", *[k for _, _, k in REFERENCE]),))
    con.commit(); con.close()
    return str(d)


# -- frame decoding ----------------------------------------------------------

def test_encode_decode_round_trip():
    rnd = np.random.default_rng(1)
    n = 3000
    scan = np.sort(rnd.integers(0, 64, n))
    tof = np.empty(n, np.int64)
    for s in np.unique(scan):
        m = scan == s
        tof[m] = np.sort(rnd.choice(400000, m.sum(), replace=False))
    inten = rnd.integers(0, 2 ** 32, n, dtype=np.uint64).astype(np.uint32)
    blk = encode_frame_payload(64, scan, tof, inten)
    import zstandard
    fr = decode_frame_payload(zstandard.ZstdDecompressor().decompress(blk[8:]))
    assert fr.num_scans == 64
    np.testing.assert_array_equal(fr.scan, scan)
    np.testing.assert_array_equal(fr.tof_indices, tof)
    np.testing.assert_array_equal(fr.intensities, inten)


def test_decoder_rejects_garbage():
    with pytest.raises(TdfBinError):
        decode_frame_payload(b"\x01\x02\x03")
    with pytest.raises(TdfBinError):
        decode_frame_payload(bytes(4 * 3))   # scan_count 0


@needs_oracle
def test_frames_match_vendor(dataset):
    with TdfData(dataset) as tdf, TdfOracle() as oracle:
        oracle.open(dataset)
        total = 0
        for fr in tdf.frames:
            ours = tdf.read_frame(fr.id)
            vs, vi, vint = oracle.read_scans(fr.id, 0, fr.num_scans)
            np.testing.assert_array_equal(ours.scan, vs)
            np.testing.assert_array_equal(ours.tof_indices, vi)
            np.testing.assert_array_equal(ours.intensities, vint)
            total += ours.num_peaks
        assert total > 4000


def test_mz_axis_recovers_the_published_equation(dataset):
    """The synthetic file was written with the OpenMS/Gobom equation at known
    coefficients; the fitted calibration must reproduce it inside the
    calibrant range. (The vendor library is not consulted for m/z: the
    synthetic MzCalibration row holds inert placeholders.)"""
    cal = synth.Calibration()
    with TdfData(dataset) as tdf:
        mz = np.linspace(tdf.calibration.mz_min, tdf.calibration.mz_max, 501)
        idx = cal.mz_to_index(mz)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ours = tdf.index_to_mz(idx)
        ppm = np.abs(ours - mz) / mz * 1e6
        assert ppm.max() < 0.01


# -- mobility --------------------------------------------------------------------

def test_reference_ions_are_found_and_calibration_is_linear(dataset):
    with warnings.catch_warnings():
        warnings.simplefilter("error")   # no fallback warning expected
        with TdfData(dataset) as tdf:
            assert tdf.mobility_fallback_reason is None
            assert tdf.mobility.method == "reference_ions"
            assert len(tdf.mobility.reference_ions) == 5
            for ion in tdf.mobility.reference_ions:
                assert ion.found and ion.snr >= 50 and abs(ion.mz_error_ppm) < 2.5
            # data were generated with 1/K0 linear in scan between the recorded range
            lo, hi = tdf.meta.oneoverk0_acq_range
            scans = np.arange(NUM_SCANS)
            truth = hi + (lo - hi) * scans / (NUM_SCANS - 1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")   # extrapolation notice, expected here
                got = tdf.scan_to_oneoverk0(scans)
            # apexes are integer scans, so expect sub-scan error
            assert np.abs(got - truth).max() < (hi - lo) / NUM_SCANS * 1.5


def test_blank_falls_back_to_acquisition_range(blank):
    with pytest.warns(RuntimeWarning, match="Falling back"):
        tdf = TdfData(blank)
    try:
        assert tdf.mobility.method == "acquisition_range"
        assert "quality gates" in tdf.mobility_fallback_reason
        lo, hi = tdf.meta.oneoverk0_acq_range
        assert tdf.scan_to_oneoverk0(0) == pytest.approx(hi)
        assert tdf.scan_to_oneoverk0(NUM_SCANS - 1) == pytest.approx(lo)
    finally:
        tdf.close()


def test_reference_ions_mode_refuses_blank(blank):
    with pytest.raises(TsfError, match="reference ions"):
        TdfData(blank, mobility="reference_ions")


def test_explicit_calibration_is_honoured(dataset):
    with pytest.warns(RuntimeWarning, match="nominal acquisition range"):
        cal = MobilityCalibration.from_acquisition_range(0.5, 1.5, NUM_SCANS)
    with TdfData(dataset, mobility=cal) as tdf:
        assert tdf.mobility is cal
        assert tdf.scan_to_oneoverk0(0) == pytest.approx(1.5)


def test_mobility_round_trip(dataset):
    with TdfData(dataset) as tdf:
        k0 = np.linspace(0.75, 1.5, 50)
        back = tdf.scan_to_oneoverk0(tdf.oneoverk0_to_scan(k0))
        np.testing.assert_allclose(back, k0, rtol=1e-9)


@pytest.mark.filterwarnings("ignore:.*outside the:RuntimeWarning")
def test_frame_mz_mobility(dataset):
    with TdfData(dataset) as tdf:
        mz, k0, inten = tdf.frame_mz_mobility(3)
        assert mz.shape == k0.shape == inten.shape
        assert np.all(np.diff(tdf.read_frame(3).scan) >= 0)


def test_pasef_tables(dataset):
    with TdfData(dataset) as tdf:
        rows = tdf.meta.pasef_msms_info(5)
        assert len(rows) == 1 and rows[0]["IsolationMz"] == pytest.approx(622.03)
        assert tdf.meta.frame(5).msms_type == 8


# -- CCS -----------------------------------------------------------------------

def test_ccs_round_trip():
    k0 = np.array([0.8, 1.0, 1.4]); mz = np.array([400.0, 622.0, 1200.0])
    ccs = oneoverk0_to_ccs(k0, mz, charge=1)
    np.testing.assert_allclose(ccs_to_oneoverk0(ccs, mz, charge=1), k0, rtol=1e-12)
    # Mason-Schamp: CCS grows with 1/K0 and with charge
    assert np.all(np.diff(ccs) > 0)
    assert oneoverk0_to_ccs(1.0, 622.0, charge=2) > oneoverk0_to_ccs(1.0, 622.0, charge=1)


@needs_oracle
def test_ccs_matches_vendor():
    for k0, z, mz in [(0.7363, 1, 322.0481), (0.9915, 1, 622.029), (1.3934, 1, 1221.9906), (1.1, 2, 900.0)]:
        with TdfOracle() as oracle:
            assert oneoverk0_to_ccs(k0, mz, z) == pytest.approx(oracle.oneoverk0_to_ccs(k0, z, mz), rel=1e-6)
