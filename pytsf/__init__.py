"""pytsf — read Bruker timsTOF data (TSF and TDF) without the Bruker SDK.

Bruker's tdf-sdk ships `timsdata.dll` (Windows) and `libtimsdata.so`
(Linux x86-64) only, so `.d` data cannot be read on macOS — including Apple
Silicon — using the vendor library. This package reimplements the reading path
in pure Python so it runs anywhere NumPy does.

    from pytsf import TsfData, TdfData

    with TsfData("/path/to/imaging.d") as tsf:          # TSF: single spectra
        mz, intensity = tsf.line_spectrum_mz(1)

    with TdfData("/path/to/pasef.d") as tdf:            # TDF: ion mobility
        mz, oneoverk0, intensity = tdf.frame_mz_mobility(1)

See README.md for how m/z and mobility are calibrated and how accurate they are.
"""

from pytsf.api import TsfData
from pytsf.meta import Frame, MzCalibrationRow, TsfError, TsfMetadata
from pytsf.mobility import MobilityCalibration, ccs_to_oneoverk0, find_reference_ions, oneoverk0_to_ccs
from pytsf.tdf import TdfData, TdfFrame, TdfMetadata

__all__ = [
    "Frame",
    "MobilityCalibration",
    "MzCalibrationRow",
    "TdfData",
    "TdfFrame",
    "TdfMetadata",
    "TsfData",
    "TsfError",
    "TsfMetadata",
    "__version__",
    "ccs_to_oneoverk0",
    "find_reference_ions",
    "oneoverk0_to_ccs",
]

__version__ = "1.1.0"
