"""Frame ids given as numpy integers must behave like plain ints.

sqlite3 binds ``np.int64`` as a BLOB, so a raw ``WHERE Id=?`` lookup with one
matches nothing. Readers that pull frame ids out of numpy arrays (e.g.
``MaldiFrameInfo.Frame``) rely on pytsf accepting them, as the vendor API does.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from pytsf import TdfData, TsfData  # noqa: E402

import synth  # noqa: E402
from test_against_vendor import _spectra  # noqa: E402
from test_tdf import NUM_SCANS, REFERENCE, _frames  # noqa: E402


@pytest.fixture
def tsf_dataset(tmp_path):
    return synth.build(str(tmp_path / "analysis.d"), _spectra())


@pytest.fixture
def tdf_dataset(tmp_path):
    return synth.build_tdf(str(tmp_path / "analysis.d"), _frames(), num_scans=NUM_SCANS, reference_ions=REFERENCE)


def test_tsf_numpy_frame_id(tsf_dataset):
    """Each TSF reader returns the same data for ``np.int64(k)`` and ``k``."""
    with TsfData(tsf_dataset) as tsf:
        ids = tsf.meta.frame_ids()
        # numpy ids first, so a plain-int cache entry cannot mask the lookup
        np_results = [
            (
                tsf.meta.frame(np.int64(k)),
                tsf.read_line_spectrum(np.int64(k)),
                tsf.read_line_spectrum_with_width(np.int64(k)),
                tsf.read_profile_spectrum(np.int64(k)),
            )
            for k in ids
        ]
    with TsfData(tsf_dataset) as tsf:
        for k, (frame, line, line_w, prof) in zip(ids, np_results):
            assert frame == tsf.meta.frame(k)
            assert type(frame.id) is int
            for a, b in zip(line, tsf.read_line_spectrum(k)):
                np.testing.assert_array_equal(a, b)
            for a, b in zip(line_w, tsf.read_line_spectrum_with_width(k)):
                np.testing.assert_array_equal(a, b)
            np.testing.assert_array_equal(prof, tsf.read_profile_spectrum(k))
            assert tsf.meta.properties(np.int64(k)) == tsf.meta.properties(k)


def test_tdf_numpy_frame_id(tdf_dataset):
    """Each TDF reader returns the same data for ``np.int64(k)`` and ``k``."""
    with TdfData(tdf_dataset) as tdf:
        ids = tdf.meta.frame_ids()
        np_results = [
            (
                tdf.meta.frame(np.int64(k)),
                tdf.read_frame(np.int64(k)),
                [tuple(r) for r in tdf.meta.pasef_msms_info(np.int64(k))],
            )
            for k in ids
        ]
    with TdfData(tdf_dataset) as tdf:
        for k, (frame, fr, pasef) in zip(ids, np_results):
            assert frame == tdf.meta.frame(k)
            ref = tdf.read_frame(k)
            np.testing.assert_array_equal(fr.scan, ref.scan)
            np.testing.assert_array_equal(fr.tof_indices, ref.tof_indices)
            np.testing.assert_array_equal(fr.intensities, ref.intensities)
            assert pasef == [tuple(r) for r in tdf.meta.pasef_msms_info(k)]
        assert len(tdf.meta.pasef_msms_info(np.int64(5))) == 1
