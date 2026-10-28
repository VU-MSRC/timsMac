"""SQLite metadata layer for Bruker TSF (timsTOF Spectrum Format) analyses.

This module reads `analysis.tsf`, which is a plain SQLite database whose schema
Bruker documents in the `tsf-schema.sql` reference shipped with the SDK. Every
table and column used here is one that reference describes.
"""

from __future__ import annotations

import os
import sqlite3
from dataclasses import dataclass
from typing import Any

TSF_SQLITE_NAME = "analysis.tsf"
TSF_BINARY_NAME = "analysis.tsf_bin"

# Compression type recorded in GlobalMetadata for TSF analyses. TDF (the
# ion-mobility sibling format) uses 2; TSF uses 3. The Bruker library hard
# refuses anything else for TSF, so we mirror that.
TSF_COMPRESSION_TYPE = 3


class TsfError(Exception):
    """Raised for malformed or unsupported TSF analyses."""


@dataclass(frozen=True)
class MzCalibrationRow:
    """The parts of a `MzCalibration` row this reader uses.

    Only the digitiser timing is needed: flight time is affine in sample
    number, `t = DigitizerDelay + index * DigitizerTimebase`, and the
    calibration itself is fitted from the reference points the file records
    (see `calibration_fitted`).
    """

    id: int
    digitizer_timebase: float
    digitizer_delay: float


@dataclass(frozen=True)
class Frame:
    """One spectrum. Bruker calls these "frames" in the schema so that the same
    queries work against TDF and TSF files alike; in TSF each frame is a single
    mass spectrum with no ion-mobility dimension.
    """

    id: int
    time: float
    polarity: str
    scan_mode: int
    msms_type: int
    tims_id: int | None
    max_intensity: int
    summed_intensities: int
    num_peaks: int | None
    mz_calibration: int
    t1: float
    t2: float
    property_group: int | None = None

    @property
    def has_data(self) -> bool:
        return self.tims_id is not None


class TsfMetadata:
    """Read-only access to the `analysis.tsf` SQLite database."""

    def __init__(self, analysis_directory: str | os.PathLike):
        self.analysis_directory = str(analysis_directory)
        if not os.path.isdir(self.analysis_directory):
            raise TsfError(f"Specified analysis directory does not exist: {self.analysis_directory}")

        self.sqlite_path = os.path.join(self.analysis_directory, TSF_SQLITE_NAME)
        self.binary_path = os.path.join(self.analysis_directory, TSF_BINARY_NAME)
        if not os.path.isfile(self.sqlite_path):
            raise TsfError(f"No {TSF_SQLITE_NAME} in {self.analysis_directory}")

        # Open read-only so we can never damage an acquisition directory.
        uri = "file:" + self.sqlite_path.replace("?", "%3f").replace("#", "%23") + "?mode=ro"
        self.conn = sqlite3.connect(uri, uri=True)
        self.conn.row_factory = sqlite3.Row

        self.global_metadata: dict[str, str] = {
            row["Key"]: row["Value"] for row in self.conn.execute("SELECT Key, Value FROM GlobalMetadata")
        }
        self._validate()

        self._frame_cache: dict[int, Frame] = {}
        self._frames_have_property_group = "PropertyGroup" in self._columns("Frames")
        self._cal_cache: dict[int, MzCalibrationRow] = {}

    # -- validation ------------------------------------------------------

    def _validate(self) -> None:
        schema_type = self.global_metadata.get("SchemaType")
        if schema_type is not None and schema_type.upper() != "TSF":
            raise TsfError(
                f"Not a TSF analysis: SchemaType={schema_type!r}. "
                "This reader does not handle TDF (ion-mobility) data."
            )

        raw = self.global_metadata.get("TimsCompressionType")
        if raw is None:
            raise TsfError("GlobalMetadata is missing TimsCompressionType.")
        try:
            self.compression_type = int(raw)
        except ValueError as exc:
            raise TsfError(f"Non-integer TimsCompressionType={raw!r}") from exc
        if self.compression_type != TSF_COMPRESSION_TYPE:
            raise TsfError(
                f"TimsCompressionType {self.compression_type} not supported " f"(TSF expects {TSF_COMPRESSION_TYPE})."
            )

    def _table_exists(self, name: str) -> bool:
        row = self.conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?",
            (name,),
        ).fetchone()
        return row is not None

    def _columns(self, table: str) -> list[str]:
        return [r[1] for r in self.conn.execute(f'PRAGMA table_info("{table}")')]

    # -- convenience metadata --------------------------------------------

    def _meta_flag(self, key: str) -> bool:
        return str(self.global_metadata.get(key, "0")).strip() in ("1", "true", "True")

    @property
    def has_line_spectra(self) -> bool:
        return self._meta_flag("HasLineSpectra")

    @property
    def has_line_spectra_peak_width(self) -> bool:
        return self._meta_flag("HasLineSpectraPeakWidth")

    @property
    def has_profile_spectra(self) -> bool:
        return self._meta_flag("HasProfileSpectra")

    @property
    def digitizer_num_samples(self) -> int | None:
        raw = self.global_metadata.get("DigitizerNumSamples")
        return int(raw) if raw is not None else None

    @property
    def mz_acq_range(self) -> tuple[float | None, float | None]:
        lo = self.global_metadata.get("MzAcqRangeLower")
        hi = self.global_metadata.get("MzAcqRangeUpper")
        return (float(lo) if lo else None, float(hi) if hi else None)

    @property
    def frame_count(self) -> int:
        return self.conn.execute("SELECT COUNT(*) FROM Frames").fetchone()[0]

    @property
    def is_maldi(self) -> bool:
        """MALDI imaging analyses carry a MaldiFrameInfo table (written by the
        acquisition software; not part of the core documented schema).
        """
        return self._table_exists("MaldiFrameInfo")

    # -- frames ----------------------------------------------------------

    def frame_ids(self) -> list[int]:
        return [r[0] for r in self.conn.execute("SELECT Id FROM Frames ORDER BY Id")]

    def frame(self, frame_id: int) -> Frame:
        cached = self._frame_cache.get(frame_id)
        if cached is not None:
            return cached
        row = self.conn.execute(
            "SELECT Id, Time, Polarity, ScanMode, MsMsType, TimsId, MaxIntensity, "
            "SummedIntensities, NumPeaks, MzCalibration, T1, T2 FROM Frames WHERE Id=?",
            (frame_id,),
        ).fetchone()
        if row is None:
            raise TsfError(f"No such spectrum id: {frame_id}")
        pg = None
        if self._frames_have_property_group:
            pg = self.conn.execute("SELECT PropertyGroup FROM Frames WHERE Id=?", (frame_id,)).fetchone()[0]
        frame = Frame(
            id=row["Id"],
            time=row["Time"],
            polarity=row["Polarity"],
            scan_mode=row["ScanMode"],
            msms_type=row["MsMsType"],
            tims_id=row["TimsId"],
            max_intensity=row["MaxIntensity"],
            summed_intensities=row["SummedIntensities"],
            num_peaks=row["NumPeaks"],
            mz_calibration=row["MzCalibration"],
            t1=row["T1"],
            t2=row["T2"],
            property_group=pg,
        )
        self._frame_cache[frame_id] = frame
        return frame

    def frames(self) -> list[Frame]:
        return [self.frame(i) for i in self.frame_ids()]

    # -- calibration -----------------------------------------------------

    def mz_calibration(self, frame_id: int) -> MzCalibrationRow:
        """The digitiser timing for one spectrum's calibration row."""
        frame = self.frame(frame_id)
        cached = self._cal_cache.get(frame.mz_calibration)
        if cached is not None:
            return cached
        row = self.conn.execute(
            "SELECT Id, DigitizerTimebase, DigitizerDelay FROM MzCalibration " "WHERE Id=?",
            (frame.mz_calibration,),
        ).fetchone()
        if row is None:
            raise TsfError(f"Spectrum {frame_id} references missing MzCalibration id " f"{frame.mz_calibration}")
        cal = MzCalibrationRow(
            id=row["Id"],
            digitizer_timebase=float(row["DigitizerTimebase"]),
            digitizer_delay=float(row["DigitizerDelay"]),
        )
        self._cal_cache[frame.mz_calibration] = cal
        return cal

    # -- generic helpers --------------------------------------------------

    def table(self, name: str) -> list[sqlite3.Row]:
        """Read any table present in this analysis (e.g. `MaldiFrameInfo`,
        `FrameMsMsInfo`, `Segments`), so files carrying vendor extensions beyond
        the documented schema remain usable.
        """
        if not self._table_exists(name):
            raise TsfError(f"No table or view named {name!r} in this analysis.")
        return list(self.conn.execute(f'SELECT * FROM "{name}"'))

    def properties(self, frame_id: int) -> dict[str, Any]:
        """Generic frame properties resolved through the documented `Properties`
        view (frame-level values overriding group-level ones).
        """
        if not self._table_exists("Properties"):
            return {}
        rows = self.conn.execute(
            "SELECT pd.PermanentName AS name, p.Value AS value "
            "FROM Properties p JOIN PropertyDefinitions pd ON p.Property = pd.Id "
            "WHERE p.Frame = ?",
            (frame_id,),
        )
        return {r["name"]: r["value"] for r in rows if r["value"] is not None}

    def close(self) -> None:
        if getattr(self, "conn", None) is not None:
            self.conn.close()
            self.conn = None  # type: ignore[assignment]

    def __enter__(self) -> TsfMetadata:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
