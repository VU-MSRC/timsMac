"""Public API for reading Bruker TSF analyses without the Bruker library.

`TsfData` deliberately mirrors the method names of Bruker's own
`examples/py/tsfdata.py` wrapper, so existing scripts can usually be ported by
changing only the import:

    # before
    import tsfdata
    tsf = tsfdata.TsfData(analysis_dir)

    # after
    from pytsf import TsfData
    tsf = TsfData(analysis_dir)

Pythonic aliases (`read_line_spectrum`, `index_to_mz`, ...) are provided
alongside the camelCase originals.
"""

from __future__ import annotations

import contextlib
import os
from collections.abc import Iterator

import numpy as np

from pytsf import bin as _bin
from pytsf import calibration_fitted as _fitted
from pytsf.meta import Frame, TsfError, TsfMetadata

__all__ = ["Frame", "TsfData", "TsfError"]


class TsfData:
    """A Bruker timsTOF .d analysis directory in TSF format.

    Parameters
    ----------
    analysis_directory
        The `.d` directory containing `analysis.tsf` and `analysis.tsf_bin`.
    use_recalibrated_state
        Accepted for signature compatibility with the Bruker wrapper. Post-hoc
        post-hoc recalibration states are not interpreted by this reader;
        passing True raises rather than silently returning values the caller
        would believe were recalibrated.
    """

    def __init__(
        self,
        analysis_directory: str | os.PathLike,
        use_recalibrated_state: bool = False,
    ) -> None:
        if use_recalibrated_state:
            raise TsfError(
                "use_recalibrated_state=True is not supported: this reader only "
                "implements the raw acquisition-time calibration. Pass False to "
                "read raw-calibration m/z values."
            )
        self.meta = TsfMetadata(analysis_directory)
        self.analysis_directory = self.meta.analysis_directory

        if not os.path.isfile(self.meta.binary_path):
            raise TsfError(f"Missing binary file: {self.meta.binary_path}")
        self._bin = _bin.BinaryFile(self.meta.binary_path)
        self._fitted = _fitted.FittedCalibration(self.meta.sqlite_path)

    @property
    def calibration(self) -> _fitted.FittedCalibration:
        """The fitted calibration. One fit covers the whole analysis."""
        return self._fitted

    # -- compatibility shim ----------------------------------------------

    @property
    def conn(self):
        """The underlying sqlite3 connection, as Bruker's wrapper exposes."""
        return self.meta.conn

    # -- m/z conversion ---------------------------------------------------

    def index_to_mz(self, frame_id: int, indices) -> np.ndarray:
        """Convert (possibly non-integer) index values to m/z for one spectrum."""
        return self._fitted.index_to_mz(indices)

    def mz_to_index(self, frame_id: int, mzs) -> np.ndarray:
        """Convert m/z values to index values for one spectrum."""
        return self._fitted.mz_to_index(mzs)

    indexToMz = index_to_mz
    mzToIndex = mz_to_index

    # -- spectra ----------------------------------------------------------

    def read_line_spectrum(self, frame_id: int):
        """Centroided ("line") spectrum as (index, intensity).

        Use `index_to_mz` to turn the indices into m/z values.
        """
        if not self.meta.has_line_spectra:
            raise TsfError("This analysis contains no line spectra " "(GlobalMetadata.HasLineSpectra == 0).")
        frame = self.meta.frame(frame_id)
        idx, inten, _ = self._bin.read_line_spectrum(frame.tims_id, frame.num_peaks)
        return idx, inten

    def read_line_spectrum_with_width(self, frame_id: int):
        """Centroided spectrum as (index, intensity, width)."""
        if not self.meta.has_line_spectra:
            raise TsfError("This analysis contains no line spectra " "(GlobalMetadata.HasLineSpectra == 0).")
        if not self.meta.has_line_spectra_peak_width:
            raise TsfError("This analysis contains no peak widths " "(GlobalMetadata.HasLineSpectraPeakWidth == 0).")
        frame = self.meta.frame(frame_id)
        return self._bin.read_line_spectrum(frame.tims_id, frame.num_peaks, want_width=True)

    def read_profile_spectrum(self, frame_id: int) -> np.ndarray:
        """Profile spectrum as a dense uint32 intensity array; the position in
        the array is the index (so `index_to_mz(frame_id, arange(len(a)))`
        gives the m/z axis).
        """
        if not self.meta.has_profile_spectra:
            raise TsfError("This analysis contains no profile spectra " "(GlobalMetadata.HasProfileSpectra == 0).")
        frame = self.meta.frame(frame_id)
        return self._bin.read_profile_spectrum(frame.tims_id)

    readLineSpectrum = read_line_spectrum
    readLineSpectrumWithWidth = read_line_spectrum_with_width
    readProfileSpectrum = read_profile_spectrum

    # -- conveniences the Bruker wrapper lacks ----------------------------

    def line_spectrum_mz(self, frame_id: int):
        """(m/z, intensity) for a centroided spectrum, in one call."""
        idx, inten = self.read_line_spectrum(frame_id)
        return self.index_to_mz(frame_id, idx), inten

    def profile_spectrum_mz(self, frame_id: int):
        """(m/z, intensity) for a profile spectrum, in one call."""
        inten = self.read_profile_spectrum(frame_id)
        # One fit per analysis, so the axis depends only on the length and is
        # cached; recomputing it per spectrum would dominate an imaging run.
        return self._fitted.index_axis(inten.size), inten

    @property
    def frames(self) -> list[Frame]:
        return self.meta.frames()

    def __len__(self) -> int:
        return self.meta.frame_count

    def __iter__(self) -> Iterator[Frame]:
        return iter(self.meta.frames())

    # -- lifecycle --------------------------------------------------------

    def close(self) -> None:
        if getattr(self, "_bin", None) is not None:
            self._bin.close()
            self._bin = None  # type: ignore[assignment]
        if getattr(self, "meta", None) is not None:
            self.meta.close()

    def __enter__(self) -> TsfData:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.close()
