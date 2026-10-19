"""pytsf — read Bruker timsTOF TSF data without the Bruker SDK.

Bruker's tdf-sdk ships `timsdata.dll` (Windows) and `libtimsdata.so`
(Linux x86-64) only, so `analysis.tsf` data cannot be read on macOS — including
Apple Silicon — using the vendor library. This package reimplements the TSF
reading path in pure Python so it runs anywhere NumPy does.

    from pytsf import TsfData

    with TsfData("/path/to/analysis.d") as tsf:
        print(len(tsf), "spectra")
        mz, intensity = tsf.line_spectrum_mz(1)

Scope: TSF only (single mass spectra, e.g. MALDI imaging and non-ion-mobility
timsTOF acquisitions). It does not read TDF, the ion-mobility sibling format --
for that, use OpenTIMS or timsrust, which handle TDF well.
"""

from .api import TsfData
from .meta import Frame, MzCalibrationRow, TsfError, TsfMetadata

__all__ = [
    "TsfData",
    "TsfMetadata",
    "Frame",
    "MzCalibrationRow",
    "TsfError",
    "__version__",
]

__version__ = "1.0.0"
