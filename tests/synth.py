"""Build synthetic TSF analyses for testing.

Writes a directory containing `analysis.tsf` (SQLite, per the schema Bruker
documents in tsf-schema.sql) and `analysis.tsf_bin` (encoded with
`pytsf.bin.encode_spectrum_block`). Used both to test `pytsf` against itself
and to feed identical input to the vendor library for differential testing.
"""

from __future__ import annotations

import os
import shutil
import sqlite3
from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

import struct

from pytsf.bin import encode_spectrum_block

SCALE = 1e12


def published_mz_to_time(ml1: float, ml2: float, ml3: float, m):
    """The published three-point TOF equation (OpenMS TOFCalibration; Gobom et
    al., Anal. Chem. 2002, 74(15), 3915-3923). Used here only to synthesise
    self-consistent test data."""
    m = np.asarray(m, dtype=float)
    return ml2 + np.sqrt(SCALE / ml1 * m) + ml3 * m

SCHEMA = """
CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
-- pytsf reads only Id / DigitizerTimebase / DigitizerDelay from this table.
-- The remaining documented columns are declared, and filled with the constant
-- 1.0, for one reason only: the vendor reference implementation refuses to
-- open a file in which they are absent, and the test suite compares spectrum
-- decoding against it. No test asks the vendor for m/z on a synthetic file.
CREATE TABLE MzCalibration (
    Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL,
    DigitizerTimebase REAL NOT NULL, DigitizerDelay REAL NOT NULL,
    T1 REAL NOT NULL, T2 REAL NOT NULL, dC1 REAL NOT NULL, dC2 REAL NOT NULL,
    C0, C1, C2, C3, C4, C5, C6, C7);
CREATE TABLE Frames (
    Id INTEGER PRIMARY KEY, Time REAL NOT NULL, Polarity CHAR(1) NOT NULL,
    ScanMode INTEGER NOT NULL, MsMsType INTEGER NOT NULL, TimsId INTEGER,
    MaxIntensity INTEGER NOT NULL, SummedIntensities INTEGER NOT NULL,
    NumPeaks INTEGER, MzCalibration INTEGER NOT NULL,
    T1 REAL NOT NULL, T2 REAL NOT NULL, PropertyGroup INTEGER);
CREATE TABLE CalibrationInfo (
    KeyPolarity CHAR(1), KeyName TEXT, Value TEXT,
    PRIMARY KEY (KeyPolarity, KeyName));
"""


@dataclass
class Calibration:
    """A synthetic calibration expressed in the published equation's own terms."""
    ml1: float = 1.5e6
    ml2: float = 1234.5
    ml3: float = 0.0
    digitizer_timebase: float = 0.08
    digitizer_delay: float = -500.0

    def mz_to_index(self, m):
        t = published_mz_to_time(self.ml1, self.ml2, self.ml3, m)
        return (t - self.digitizer_delay) / self.digitizer_timebase


@dataclass
class Spectrum:
    index: Sequence[float] = field(default_factory=list)
    intensity: Sequence[float] = field(default_factory=list)
    width: Optional[Sequence[float]] = None
    profile: Optional[Sequence[int]] = None
    t1: float = 25.0
    t2: float = 30.0
    gap_before: int = 0


def build(
    outdir: str,
    spectra: Sequence[Spectrum],
    calibration: Calibration | None = None,
    num_samples: int = 400000,
    has_line: bool = True,
    has_profile: bool = True,
    has_width: bool | None = None,
    acquisition_software: str = "synthetic",
    reference_points: bool = True,
) -> str:
    """Write a synthetic TSF analysis directory and return its path."""
    cal = calibration or Calibration()
    if has_width is None:
        has_width = any(s.width is not None for s in spectra)

    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    blob = bytearray()
    offsets: list[Optional[int]] = []
    for s in spectra:
        blob += b"\x00" * s.gap_before
        if not s.index and s.profile is None:
            offsets.append(None)   # TimsId NULL: "spectrum has no data"
            continue
        offsets.append(len(blob))
        blob += encode_spectrum_block(
            s.index, s.intensity, s.width, s.profile
        )

    with open(os.path.join(outdir, "analysis.tsf_bin"), "wb") as f:
        f.write(bytes(blob))

    con = sqlite3.connect(os.path.join(outdir, "analysis.tsf"))
    con.executescript(SCHEMA)
    con.executemany(
        "INSERT INTO GlobalMetadata VALUES (?,?)",
        {
            "SchemaType": "TSF",
            "SchemaVersionMajor": "3",
            "SchemaVersionMinor": "4",
            # TSF is compression type 3; the vendor library refuses any other
            # value for TSF.
            "TimsCompressionType": "3",
            "DigitizerNumSamples": str(num_samples),
            "HasLineSpectra": "1" if has_line else "0",
            "HasProfileSpectra": "1" if has_profile else "0",
            "HasLineSpectraPeakWidth": "1" if has_width else "0",
            "MzAcqRangeLower": "100",
            "MzAcqRangeUpper": "1700",
            "AcquisitionSoftware": acquisition_software,
            "AcquisitionSoftwareVersion": "0",
            "InstrumentName": "synthetic",
        }.items(),
    )
    con.execute(
        "INSERT INTO MzCalibration (Id,ModelType,DigitizerTimebase,"
        "DigitizerDelay,T1,T2,dC1,dC2,C0,C1,C2,C3,C4,C5,C6,C7) "
        "VALUES (1,1,?,?,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)",
        (cal.digitizer_timebase, cal.digitizer_delay),
    )
    for i, (s, off) in enumerate(zip(spectra, offsets), start=1):
        inten = np.asarray(s.intensity, dtype=np.float64)
        con.execute(
            "INSERT INTO Frames (Id,Time,Polarity,ScanMode,MsMsType,TimsId,"
            "MaxIntensity,SummedIntensities,NumPeaks,MzCalibration,T1,T2,"
            "PropertyGroup) VALUES (?,?,'+',0,0,?,?,?,?,1,?,?,NULL)",
            (i, i * 0.1, off,
             int(inten.max()) if inten.size else 0,
             int(inten.sum()) if inten.size else 0,
             len(s.index), s.t1, s.t2),
        )
    if reference_points:
        # Reference calibration points, as a real acquisition records them,
        # generated here from the published equation so the fitted calibration has
        # something self-consistent to fit.
        ref_mz = np.array([400.0, 600.0, 850.0, 1100.0, 1400.0, 1700.0])
        ref_t = published_mz_to_time(cal.ml1, cal.ml2, cal.ml3, ref_mz)
        con.executemany("INSERT INTO CalibrationInfo VALUES (?,?,?)", [
            ("+", "ReferenceMassList", "synthetic"),
            ("+", "ReferencePeakMasses",
             struct.pack(f"<{ref_mz.size}d", *ref_mz)),
            ("+", "MeasuredTimesOfFlight",
             struct.pack(f"<{ref_t.size}d", *ref_t)),
        ])

    con.commit()
    con.close()
    return outdir


# ---------------------------------------------------------------------------
# TDF (ion mobility) synthetic analyses
# ---------------------------------------------------------------------------

from pytsf.tdf_bin import encode_frame_payload  # noqa: E402

TDF_SCHEMA = """
CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
CREATE TABLE MzCalibration (
    Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL,
    DigitizerTimebase REAL NOT NULL, DigitizerDelay REAL NOT NULL,
    T1 REAL NOT NULL, T2 REAL NOT NULL, dC1 REAL NOT NULL, dC2 REAL NOT NULL,
    C0, C1, C2, C3, C4, C5, C6, C7);
CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL,
    C0, C1, C2, C3, C4, C5, C6, C7, C8, C9);
CREATE TABLE Frames (
    Id INTEGER PRIMARY KEY, Time REAL NOT NULL, Polarity CHAR(1) NOT NULL,
    ScanMode INTEGER NOT NULL, MsMsType INTEGER NOT NULL, TimsId INTEGER,
    MaxIntensity INTEGER NOT NULL, SummedIntensities INTEGER NOT NULL,
    NumScans INTEGER NOT NULL, NumPeaks INTEGER NOT NULL,
    MzCalibration INTEGER NOT NULL, T1 REAL NOT NULL, T2 REAL NOT NULL,
    TimsCalibration INTEGER NOT NULL, PropertyGroup INTEGER,
    AccumulationTime REAL NOT NULL, RampTime REAL NOT NULL, Pressure REAL);
CREATE TABLE PasefFrameMsMsInfo (
    Frame INTEGER NOT NULL, ScanNumBegin INTEGER NOT NULL, ScanNumEnd INTEGER NOT NULL,
    IsolationMz REAL NOT NULL, IsolationWidth REAL NOT NULL, CollisionEnergy REAL NOT NULL,
    Precursor INTEGER, PRIMARY KEY (Frame, ScanNumBegin));
CREATE TABLE Precursors (Id INTEGER PRIMARY KEY, LargestPeakMz REAL NOT NULL,
    AverageMz REAL NOT NULL, MonoisotopicMz REAL, Charge INTEGER, ScanNumber REAL NOT NULL,
    Intensity REAL NOT NULL, Parent INTEGER NOT NULL);
CREATE TABLE CalibrationInfo (
    KeyPolarity CHAR(1), KeyName TEXT, Value TEXT,
    PRIMARY KEY (KeyPolarity, KeyName));
"""


@dataclass
class TdfFrameSpec:
    """Peaks of one synthetic frame, as (scan, m/z, intensity) triples."""
    scan: Sequence[int]
    mz: Sequence[float]
    intensity: Sequence[int]
    msms_type: int = 0
    pasef: Optional[list] = None   # list of (scan_begin, scan_end, iso_mz, iso_width, ce, precursor)


def build_tdf(outdir: str, frames: Sequence[TdfFrameSpec], num_scans: int = 500,
              calibration: Calibration | None = None,
              k0_range: tuple = (0.6, 1.6),
              reference_ions: Optional[list] = None,
              num_samples: int = 400000) -> str:
    """Write a synthetic TDF analysis.

    The mobility axis used to *generate* the data is linear in scan between
    ``k0_range`` (high 1/K0 at scan 0). ``reference_ions`` is a list of
    (name, m/z, 1/K0); if given they are recorded in CalibrationInfo the way a
    real acquisition records its calibrant, and a peak for each is written
    into every MS1 frame at the scan its 1/K0 implies.
    """
    cal = calibration or Calibration()
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    k_hi, k_lo = k0_range[1], k0_range[0]

    def k0_to_scan(k0):
        return (k0 - k_hi) / (k_lo - k_hi) * (num_scans - 1)

    blob = bytearray()
    rows = []
    pasef_rows = []
    for i, fs in enumerate(frames, start=1):
        scan = list(fs.scan); mz = list(fs.mz); inten = list(fs.intensity)
        tof = list(np.rint(cal.mz_to_index(np.asarray(mz, dtype=float))).astype(np.int64))
        if reference_ions and fs.msms_type == 0:
            rnd = np.random.default_rng(i)
            for _, rmz, rk0 in reference_ions:
                s0 = int(round(k0_to_scan(rk0)))
                x = float(cal.mz_to_index(rmz)); x0 = int(np.floor(x)); frac = x - x0
                for ds, w in ((-1, 0.5), (0, 1.0), (1, 0.5)):
                    if 0 <= s0 + ds < num_scans:
                        # split the ion over the two neighbouring TOF bins so its
                        # intensity-weighted index is the exact fractional index
                        for xi, wi in ((x0, 1.0 - frac), (x0 + 1, frac)):
                            if wi > 0:
                                scan.append(s0 + ds); tof.append(xi); inten.append(int(round(20000 * w * wi)))
                # a sprinkle of low-level background at that m/z, as real data has
                for sb in rnd.integers(0, num_scans, 40):
                    scan.append(int(sb)); tof.append(x0 + int(rnd.integers(-2, 3))); inten.append(int(rnd.integers(1, 30)))
        tof = np.asarray(tof, dtype=np.int64)
        order = np.lexsort((tof, np.asarray(scan)))
        scan_a = np.asarray(scan)[order]; tof_a = tof[order]; int_a = np.asarray(inten, dtype=np.uint32)[order]
        # drop exact duplicates within a scan (delta would be 0 -> ambiguous)
        keep = np.ones(scan_a.size, bool)
        keep[1:] = ~((scan_a[1:] == scan_a[:-1]) & (tof_a[1:] == tof_a[:-1]))
        scan_a, tof_a, int_a = scan_a[keep], tof_a[keep], int_a[keep]
        off = len(blob)
        blob += encode_frame_payload(num_scans, scan_a, tof_a, int_a)
        rows.append((i, i * 0.1, fs.msms_type, off, int(int_a.max()) if int_a.size else 0,
                     int(int_a.sum()) if int_a.size else 0, num_scans, int(int_a.size)))
        for p in (fs.pasef or []):
            pasef_rows.append((i, *p))

    with open(os.path.join(outdir, "analysis.tdf_bin"), "wb") as f:
        f.write(bytes(blob))

    con = sqlite3.connect(os.path.join(outdir, "analysis.tdf"))
    con.executescript(TDF_SCHEMA)
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)", {
        "SchemaType": "TDF", "SchemaVersionMajor": "3", "SchemaVersionMinor": "8",
        "TimsCompressionType": "2", "DigitizerNumSamples": str(num_samples),
        "MzAcqRangeLower": "100", "MzAcqRangeUpper": "1700",
        "OneOverK0AcqRangeLower": str(k_lo), "OneOverK0AcqRangeUpper": str(k_hi),
        "AcquisitionSoftware": "synthetic", "InstrumentName": "synthetic",
        "MaxNumPeaksPerScan": "500", "PeakListIndexScaleFactor": "1",
    }.items())
    con.execute("INSERT INTO MzCalibration (Id,ModelType,DigitizerTimebase,DigitizerDelay,T1,T2,dC1,dC2,"
                "C0,C1,C2,C3,C4,C5,C6,C7) VALUES (1,1,?,?,1,1,1,1,1,1,1,1,1,1,1,1)",
                (cal.digitizer_timebase, cal.digitizer_delay))
    con.execute("INSERT INTO TimsCalibration (Id,ModelType,C0,C1,C2,C3,C4,C5,C6,C7,C8,C9) "
                "VALUES (1,1,1,1,1,1,1,1,1,1,1,1)")
    con.executemany("INSERT INTO Frames (Id,Time,Polarity,ScanMode,MsMsType,TimsId,MaxIntensity,"
                    "SummedIntensities,NumScans,NumPeaks,MzCalibration,T1,T2,TimsCalibration,"
                    "PropertyGroup,AccumulationTime,RampTime,Pressure) "
                    "VALUES (?,?,'+',8,?,?,?,?,?,?,1,25,30,1,NULL,100.0,100.0,2.4)", rows)
    con.executemany("INSERT INTO PasefFrameMsMsInfo VALUES (?,?,?,?,?,?,?)", pasef_rows)
    ref_mz = np.array([400.0, 600.0, 850.0, 1100.0, 1400.0, 1700.0])
    ref_t = published_mz_to_time(cal.ml1, cal.ml2, cal.ml3, ref_mz)
    info = [("+", "ReferencePeakMasses", struct.pack(f"<{ref_mz.size}d", *ref_mz)),
            ("+", "MeasuredTimesOfFlight", struct.pack(f"<{ref_t.size}d", *ref_t))]
    if reference_ions:
        names = b"".join(n.encode() + b"\x00" for n, _, _ in reference_ions)
        rk = np.array([k for _, _, k in reference_ions])
        # overwrite the masses with the calibrant's own masses, as a real file does
        info[0] = ("+", "ReferencePeakMasses",
                   struct.pack(f"<{len(reference_ions)}d", *[m for _, m, _ in reference_ions]))
        info[1] = ("+", "MeasuredTimesOfFlight",
                   struct.pack(f"<{len(reference_ions)}d",
                               *published_mz_to_time(cal.ml1, cal.ml2, cal.ml3,
                                                     np.array([m for _, m, _ in reference_ions]))))
        info += [("+", "ReferenceMobilityPeakNames", names),
                 ("+", "ReferencePeakMobilities", struct.pack(f"<{rk.size}d", *rk))]
    con.executemany("INSERT INTO CalibrationInfo VALUES (?,?,?)", info)
    con.commit(); con.close()
    return outdir
