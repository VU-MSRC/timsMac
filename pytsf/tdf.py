"""Reader for Bruker TDF analyses: ion-mobility (TIMS) data.

A TDF frame holds many *scans* (mobility bins); each scan holds centroided
(TOF index, intensity) peaks. This class decodes frames (``tdf_bin``), converts
TOF index to m/z with the same fitted calibration as TSF
(``calibration_fitted``), and converts scan number to 1/K0 and CCS
(``mobility``). The PASEF / DIA / PRM tables documented in Bruker's
``tdf-schema.sql`` are exposed as rows.

    from pytsf import TdfData

    with TdfData("run.d") as tdf:
        fr = tdf.read_frame(1)                 # scan, tof, intensity arrays
        mz = tdf.index_to_mz(fr.tof_indices)
        k0 = tdf.scan_to_oneoverk0(fr.scan)
        ccs = tdf.oneoverk0_to_ccs(k0, mz, charge=1)

Mobility calibration
--------------------
``mobility="auto"`` (default) looks in the data for the calibrant the file
records in ``CalibrationInfo`` (``ReferencePeakMasses`` /
``ReferencePeakMobilities``). If found, the mobility axis is fitted through the
observed ions -- the standard external calibration, ~0.3 % against the vendor.
If not found (e.g. a blank), it falls back to interpolating the documented
acquisition range and warns, since that is only good to a few percent. Pass
``mobility="acquisition_range"`` to force the fallback, or a
``MobilityCalibration`` you built yourself from your own reference ions.

Intensities
-----------
Intensities are returned exactly as stored in the file, as timsrust and
OpenTIMS also do. The vendor library returns them multiplied by a factor that
is not documented anywhere Bruker distributes; on our test files it was 1.0
when ``Frames.AccumulationTime`` was 100 ms and 2.5 when it was 40 ms. Relative
intensities within a file are unaffected. If you compare absolute intensities
with vendor-derived numbers, or across runs with different accumulation times,
be aware of this; dividing by ``Frames.AccumulationTime`` is the ordinary
physics-based normalisation and uses only a documented column.
"""

from __future__ import annotations

import os
import sqlite3
import struct
import warnings
from dataclasses import dataclass
from typing import Iterator, Optional, Union

import numpy as np

from . import tdf_bin as _bin
from .calibration_fitted import FittedCalibration
from .meta import TsfError
from .mobility import (MobilityCalibration, ccs_to_oneoverk0, find_reference_ions,
                       oneoverk0_to_ccs)

__all__ = ["TdfData", "TdfFrame", "TdfMetadata"]

TDF_SQLITE_NAME = "analysis.tdf"
TDF_BINARY_NAME = "analysis.tdf_bin"
TDF_COMPRESSION_TYPE = 2


@dataclass(frozen=True)
class TdfFrame:
    id: int
    time: float
    polarity: str
    scan_mode: int
    msms_type: int
    tims_id: Optional[int]
    num_scans: int
    num_peaks: int
    max_intensity: int
    summed_intensities: int
    accumulation_time: float
    ramp_time: float
    mz_calibration: int
    tims_calibration: int

    @property
    def has_data(self) -> bool:
        return self.tims_id is not None


class TdfMetadata:
    """Read-only access to ``analysis.tdf``."""

    def __init__(self, analysis_directory):
        self.analysis_directory = str(analysis_directory)
        if not os.path.isdir(self.analysis_directory):
            raise TsfError(f"Specified analysis directory does not exist: {self.analysis_directory}")
        self.sqlite_path = os.path.join(self.analysis_directory, TDF_SQLITE_NAME)
        self.binary_path = os.path.join(self.analysis_directory, TDF_BINARY_NAME)
        if not os.path.isfile(self.sqlite_path):
            raise TsfError(f"No {TDF_SQLITE_NAME} in {self.analysis_directory}")
        uri = "file:" + self.sqlite_path.replace("?", "%3f").replace("#", "%23") + "?mode=ro"
        self.conn = sqlite3.connect(uri, uri=True)
        self.conn.row_factory = sqlite3.Row
        self.global_metadata = {r["Key"]: r["Value"] for r in
                                self.conn.execute("SELECT Key, Value FROM GlobalMetadata")}
        st = self.global_metadata.get("SchemaType")
        if st is not None and st.upper() != "TDF":
            raise TsfError(f"Not a TDF analysis: SchemaType={st!r}. Use TsfData for TSF files.")
        ct = int(self.global_metadata.get("TimsCompressionType", "0"))
        if ct != TDF_COMPRESSION_TYPE:
            raise TsfError(
                f"TimsCompressionType {ct} not supported: this reader implements the "
                f"zstd-based type {TDF_COMPRESSION_TYPE} that current instruments write.")
        self._frames: dict[int, TdfFrame] = {}

    # -- scalar metadata --
    def _meta_float(self, key):
        v = self.global_metadata.get(key)
        return float(v) if v is not None else None

    @property
    def oneoverk0_acq_range(self):
        return (self._meta_float("OneOverK0AcqRangeLower"), self._meta_float("OneOverK0AcqRangeUpper"))

    @property
    def mz_acq_range(self):
        return (self._meta_float("MzAcqRangeLower"), self._meta_float("MzAcqRangeUpper"))

    @property
    def digitizer_num_samples(self) -> Optional[int]:
        v = self.global_metadata.get("DigitizerNumSamples")
        return int(v) if v is not None else None

    @property
    def num_scans(self) -> int:
        r = self.conn.execute("SELECT MIN(NumScans), MAX(NumScans) FROM Frames").fetchone()
        if r[0] != r[1]:
            raise TsfError(f"Frames have differing NumScans ({r[0]}..{r[1]}); not supported.")
        return int(r[0])

    @property
    def frame_count(self) -> int:
        return self.conn.execute("SELECT COUNT(*) FROM Frames").fetchone()[0]

    @property
    def is_maldi(self) -> bool:
        return self._table_exists("MaldiFrameInfo")

    def _table_exists(self, name: str) -> bool:
        return self.conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?", (name,)
        ).fetchone() is not None

    # -- frames --
    def frame_ids(self) -> list[int]:
        return [r[0] for r in self.conn.execute("SELECT Id FROM Frames ORDER BY Id")]

    def frame(self, frame_id: int) -> TdfFrame:
        f = self._frames.get(frame_id)
        if f is not None:
            return f
        r = self.conn.execute(
            "SELECT Id, Time, Polarity, ScanMode, MsMsType, TimsId, NumScans, NumPeaks, "
            "MaxIntensity, SummedIntensities, AccumulationTime, RampTime, MzCalibration, "
            "TimsCalibration FROM Frames WHERE Id=?", (frame_id,)).fetchone()
        if r is None:
            raise TsfError(f"No such frame id: {frame_id}")
        f = TdfFrame(r["Id"], r["Time"], r["Polarity"], r["ScanMode"], r["MsMsType"], r["TimsId"],
                     r["NumScans"], r["NumPeaks"], r["MaxIntensity"], r["SummedIntensities"],
                     r["AccumulationTime"], r["RampTime"], r["MzCalibration"], r["TimsCalibration"])
        self._frames[frame_id] = f
        return f

    def frames(self) -> list[TdfFrame]:
        return [self.frame(i) for i in self.frame_ids()]

    # -- documented MS/MS tables, as rows --
    def table(self, name: str) -> list[sqlite3.Row]:
        if not self._table_exists(name):
            raise TsfError(f"No table or view named {name!r} in this analysis.")
        return list(self.conn.execute(f'SELECT * FROM "{name}"'))

    def pasef_msms_info(self, frame_id: Optional[int] = None) -> list[sqlite3.Row]:
        """Rows of `PasefFrameMsMsInfo` (scan ranges, isolation, precursor id)."""
        if frame_id is None:
            return self.table("PasefFrameMsMsInfo")
        return list(self.conn.execute("SELECT * FROM PasefFrameMsMsInfo WHERE Frame=?", (frame_id,)))

    def precursors(self) -> list[sqlite3.Row]:
        return self.table("Precursors")

    def reference_calibrant(self):
        """The calibrant the file records: (names, m/z, 1/K0), or None."""
        def info(key):
            r = self.conn.execute("SELECT Value FROM CalibrationInfo WHERE KeyName=?", (key,)).fetchone()
            return r[0] if r else None
        mz, k0, names = info("ReferencePeakMasses"), info("ReferencePeakMobilities"), info("ReferenceMobilityPeakNames")
        if not (isinstance(mz, bytes) and isinstance(k0, bytes)):
            return None
        mzs = struct.unpack(f"<{len(mz) // 8}d", mz)
        k0s = struct.unpack(f"<{len(k0) // 8}d", k0)
        n = min(len(mzs), len(k0s))
        if n == 0:
            return None
        nm = [s.decode(errors="replace") for s in names.split(b"\x00") if s] if isinstance(names, bytes) else []
        nm = (nm + [f"{m:.4f}" for m in mzs])[:n]
        return nm, np.array(mzs[:n]), np.array(k0s[:n])

    def close(self):
        if getattr(self, "conn", None) is not None:
            self.conn.close(); self.conn = None  # type: ignore[assignment]


class TdfData:
    """A Bruker `.d` analysis in TDF (ion-mobility) format."""

    def __init__(self, analysis_directory, mobility: Union[str, MobilityCalibration] = "auto",
                 reference_ppm: float = 15.0, min_reference_ions: int = 3,
                 reference_quality: Optional[dict] = None):
        """
        reference_quality
            Gates a located reference ion must pass before it is used, as a
            dict overriding any of: ``max_mz_ppm`` (observed vs known m/z,
            default 2.5), ``min_snr`` (apex over the median background outside the peak,
            default 50), ``max_fit_residual_percent`` (of the final fit,
            default 1.0). Real calibrant lies on a line to ~0.15 %; a blank
            yields noise picks with ppm errors of several ppm and residuals
            of tens of percent, which these gates reject so that the reader
            falls back to the acquisition range instead of calibrating on
            noise.
        """
        self.meta = TdfMetadata(analysis_directory)
        self.analysis_directory = self.meta.analysis_directory
        if not os.path.isfile(self.meta.binary_path):
            raise TsfError(f"Missing binary file: {self.meta.binary_path}")
        self._bin = _bin.TdfBinaryFile(self.meta.binary_path)
        self._cal = FittedCalibration(self.meta.sqlite_path)
        self.reference_ions = []
        self.mobility_fallback_reason: Optional[str] = None
        q = dict(max_mz_ppm=2.5, min_snr=50.0, max_fit_residual_percent=1.0)
        q.update(reference_quality or {})
        self._mob = self._build_mobility(mobility, reference_ppm, min_reference_ions, q)

    # -- mobility calibration selection -----------------------------------

    def _build_mobility(self, mobility, ppm, min_ions, q) -> MobilityCalibration:
        ns = self.meta.num_scans
        if isinstance(mobility, MobilityCalibration):
            return mobility
        if mobility not in ("auto", "reference_ions", "acquisition_range"):
            raise TsfError("mobility must be 'auto', 'reference_ions', 'acquisition_range' "
                           "or a MobilityCalibration")
        reason = None
        if mobility in ("auto", "reference_ions"):
            cal = self.meta.reference_calibrant()
            if cal is None:
                reason = "the file records no reference calibrant in CalibrationInfo"
            else:
                names, mzs, k0s = cal
                found = find_reference_ions(self, mzs, k0s, names, ppm=ppm)
                self.reference_ions = found
                # gate 1: each ion must be a real, correctly identified peak
                ok = [f for f in found if f.found and f.snr >= q["min_snr"]
                      and abs(f.mz_error_ppm) <= q["max_mz_ppm"]]
                if len(ok) < min_ions:
                    reason = (f"only {len(ok)} of {len(found)} reference ions pass the "
                              f"quality gates (S/N >= {q['min_snr']:g}, |m/z error| <= "
                              f"{q['max_mz_ppm']:g} ppm); need {min_ions}")
                else:
                    mob = MobilityCalibration.from_reference_ions(
                        [f.apex_scan for f in ok], [f.oneoverk0 for f in ok], ns,
                        reference_ions=ok)
                    # gate 2: the fit must look like a calibration -- 1/K0 falls
                    # with scan number, and the ions lie near one line
                    slope = mob.coefficients[-2]
                    if slope >= 0:
                        reason = "reference-ion apexes do not fall with 1/K0 as TIMS requires"
                    elif mob.fit_residual_percent > q["max_fit_residual_percent"]:
                        reason = (f"reference-ion fit residual {mob.fit_residual_percent:.2f}% "
                                  f"exceeds {q['max_fit_residual_percent']:g}%")
                    else:
                        return mob
            if mobility == "reference_ions":
                raise TsfError(f"Cannot calibrate mobility from reference ions: {reason}.")
        self.mobility_fallback_reason = reason
        lo, hi = self.meta.oneoverk0_acq_range
        if lo is None or hi is None:
            raise TsfError("No OneOverK0AcqRange in GlobalMetadata; cannot calibrate mobility.")
        msg = (f"Falling back to the acquisition-range mobility approximation: {reason}."
               if reason else "Mobility calibrated from the nominal acquisition range only.")
        warnings.warn(msg + " Expect 1/K0 errors of a few percent; reference ions give ~0.3 %.",
                      RuntimeWarning, stacklevel=3)
        return MobilityCalibration.from_acquisition_range(lo, hi, ns, warn=False)

    # -- properties -------------------------------------------------------

    @property
    def calibration(self) -> FittedCalibration:
        return self._cal

    @property
    def mobility(self) -> MobilityCalibration:
        return self._mob

    @property
    def frames(self) -> list[TdfFrame]:
        return self.meta.frames()

    def __len__(self) -> int:
        return self.meta.frame_count

    def __iter__(self) -> Iterator[TdfFrame]:
        return iter(self.meta.frames())

    # -- data -------------------------------------------------------------

    def read_frame(self, frame_id: int) -> _bin.Frame:
        """Decode one frame: ``.scan``, ``.tof_indices``, ``.intensities``."""
        f = self.meta.frame(frame_id)
        fr = self._bin.read_frame(f.tims_id)
        if f.tims_id is not None and fr.num_peaks != f.num_peaks:
            raise TsfError(f"Frame {frame_id}: decoded {fr.num_peaks} peaks, Frames.NumPeaks says {f.num_peaks}")
        return fr

    def frame_mz_mobility(self, frame_id: int):
        """One frame as (m/z, 1/K0, intensity) arrays."""
        fr = self.read_frame(frame_id)
        return (self.index_to_mz(fr.tof_indices.astype(np.float64)),
                self.scan_to_oneoverk0(fr.scan.astype(np.float64)),
                fr.intensities)

    # -- conversions ------------------------------------------------------

    def index_to_mz(self, indices):
        return self._cal.index_to_mz(indices)

    def mz_to_index(self, mzs):
        return self._cal.mz_to_index(mzs)

    def scan_to_oneoverk0(self, scans):
        return self._mob.scan_to_oneoverk0(scans)

    def oneoverk0_to_scan(self, oneoverk0):
        return self._mob.oneoverk0_to_scan(oneoverk0)

    @staticmethod
    def oneoverk0_to_ccs(oneoverk0, mz, charge: int = 1):
        return oneoverk0_to_ccs(oneoverk0, mz, charge)

    @staticmethod
    def ccs_to_oneoverk0(ccs, mz, charge: int = 1):
        return ccs_to_oneoverk0(ccs, mz, charge)

    # -- lifecycle --------------------------------------------------------

    def close(self):
        if getattr(self, "_bin", None) is not None:
            self._bin.close(); self._bin = None  # type: ignore[assignment]
        if getattr(self, "meta", None) is not None:
            self.meta.close()

    def __enter__(self): return self
    def __exit__(self, *exc): self.close()
    def __del__(self):
        try: self.close()
        except Exception: pass
