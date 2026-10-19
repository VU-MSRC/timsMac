#!/usr/bin/env python3
"""Generate a small synthetic MALDI imaging TSF dataset.

The notebook uses this so that the worked example is self-contained and
reproducible without anybody's real data. The output is a valid `.d` directory
that pytsf reads; the documented `MzCalibration` columns are all present
(holding placeholder values) so that the vendor library will open it too.

    python make_demo_dataset.py [output.d]

It is synthetic: peak positions and spatial patterns are invented, and the
calibration is a plausible but made-up one. It is for exercising the code
path, not for anything scientific.
"""

from __future__ import annotations

import os
import shutil
import sqlite3
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pytsf.bin import encode_spectrum_block  # noqa: E402

SCALE = 1e12


def published_mz_to_time(ml1, ml2, ml3, m):
    """The published three-point TOF equation (OpenMS TOFCalibration; Gobom et
    al., Anal. Chem. 2002, 74(15), 3915-3923)."""
    m = np.asarray(m, dtype=float)
    return ml2 + np.sqrt(SCALE / ml1 * m) + ml3 * m

SCHEMA = """
CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
-- pytsf reads only Id / DigitizerTimebase / DigitizerDelay. The other
-- documented columns are declared and set to 1.0 so the vendor library will
-- open the file; nothing reads them.
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
CREATE TABLE MaldiFrameInfo (
    Frame INTEGER PRIMARY KEY, Chip INTEGER, SpotName TEXT,
    RegionNumber INTEGER, XIndexPos INTEGER, YIndexPos INTEGER,
    LaserPower REAL, NumLaserShots INTEGER, LaserRepRate REAL,
    MotorPositionX REAL, MotorPositionY REAL, MotorPositionZ REAL);
CREATE TABLE CalibrationInfo (
    KeyPolarity CHAR(1), KeyName TEXT, Value TEXT,
    PRIMARY KEY (KeyPolarity, KeyName));
"""

# A plausible synthetic calibration, in the published equation's own terms.
CAL = dict(timebase=0.125, delay=44000.0,
           ml1=152000.0, ml2=44000.0 + 18000.0, ml3=-4.0e-4)

WIDTH, HEIGHT = 48, 32
NUM_SAMPLES = 900000


def spatial_patterns(w: int, h: int) -> dict[float, np.ndarray]:
    """Three invented ion distributions, so the ion images show structure."""
    yy, xx = np.mgrid[0:h, 0:w]
    xn, yn = xx / (w - 1), yy / (h - 1)

    # a diagonal gradient
    grad = 0.25 + 0.75 * (xn + yn) / 2.0
    # two blobs
    blobs = (np.exp(-(((xn - 0.3) ** 2 + (yn - 0.35) ** 2) / 0.012))
             + np.exp(-(((xn - 0.72) ** 2 + (yn - 0.65) ** 2) / 0.008)))
    blobs = 0.1 + 0.9 * blobs / blobs.max()
    # a ring
    r = np.sqrt((xn - 0.5) ** 2 + (yn - 0.5) ** 2)
    ring = np.exp(-((r - 0.3) ** 2) / 0.004)
    ring = 0.1 + 0.9 * ring / ring.max()

    return {400.2345: grad, 522.3456: blobs, 744.5678: ring}


def main(argv: list[str]) -> int:
    out = argv[1] if len(argv) > 1 else "demo_imaging.d"
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)

    rng = np.random.default_rng(20260219)
    def mz_to_index(m):
        t = published_mz_to_time(CAL["ml1"], CAL["ml2"], CAL["ml3"], m)
        return (t - CAL["delay"]) / CAL["timebase"]

    patterns = spatial_patterns(WIDTH, HEIGHT)
    # a matrix background: evenly spaced peaks, spatially flat
    background_mz = np.arange(320.0, 900.0, 11.3)

    blob = bytearray()
    frames, maldi = [], []
    fid = 0
    for y in range(HEIGHT):
        for x in range(WIDTH):
            fid += 1
            mzs = list(background_mz + rng.normal(0, 0.002, background_mz.size))
            ints = list(rng.uniform(40, 160, background_mz.size))
            for mz, pat in patterns.items():
                mzs.append(mz + rng.normal(0, 0.002))
                ints.append(120.0 + 1500.0 * float(pat[y, x])
                            * rng.uniform(0.85, 1.15))
            order = np.argsort(mzs)
            mzs = np.asarray(mzs)[order]
            ints = np.round(np.asarray(ints)[order])

            idx = mz_to_index(mzs)
            offset = len(blob)
            blob += encode_spectrum_block(idx, ints, None, None)

            frames.append((fid, fid * 0.05, offset, int(ints.max()),
                           int(ints.sum()), len(idx)))
            maldi.append((fid, 0, f"R00X{x:03d}Y{y:03d}", 0, x + 1, y + 1,
                          30.0, 100, 10000.0, 20.0 * x, -20.0 * y, 0.0))

    with open(os.path.join(out, "analysis.tsf_bin"), "wb") as f:
        f.write(bytes(blob))

    con = sqlite3.connect(os.path.join(out, "analysis.tsf"))
    con.executescript(SCHEMA)
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)", {
        "SchemaType": "TSF", "SchemaVersionMajor": "3", "SchemaVersionMinor": "4",
        "TimsCompressionType": "3", "DigitizerNumSamples": str(NUM_SAMPLES),
        "HasLineSpectra": "1", "HasProfileSpectra": "0",
        "HasLineSpectraPeakWidth": "0",
        "MzAcqRangeLower": "300.000000", "MzAcqRangeUpper": "1000.000000",
        "AcquisitionSoftware": "pytsf-demo", "AcquisitionSoftwareVersion": "0",
        "InstrumentName": "synthetic", "MaldiApplicationType": "Imaging",
    }.items())
    con.execute(
        "INSERT INTO MzCalibration (Id,ModelType,DigitizerTimebase,"
        "DigitizerDelay,T1,T2,dC1,dC2,C0,C1,C2,C3,C4,C5,C6,C7) "
        "VALUES (1,1,?,?,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)",
        (CAL["timebase"], CAL["delay"]))
    con.executemany(
        "INSERT INTO Frames (Id,Time,Polarity,ScanMode,MsMsType,TimsId,"
        "MaxIntensity,SummedIntensities,NumPeaks,MzCalibration,T1,T2,"
        "PropertyGroup) VALUES (?,?,'+',20,0,?,?,?,?,1,?,?,NULL)",
        [(f[0], f[1], f[2], f[3], f[4], f[5], 25.5, 23.0) for f in frames])
    con.executemany(
        "INSERT INTO MaldiFrameInfo VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", maldi)
    # Reference calibration points, as a real acquisition records them: the
    # known calibrant masses and the flight times they were measured at. pytsf
    # fits the published TOF equation to these, so the demo exercises the same
    # path as real data. The range brackets every demo peak with margin.
    import struct
    ref_mz = np.array([310.0, 400.0, 500.0, 650.0, 800.0, 920.0])
    ref_t = published_mz_to_time(CAL["ml1"], CAL["ml2"], CAL["ml3"], ref_mz)
    con.executemany("INSERT INTO CalibrationInfo VALUES (?,?,?)", [
        ("+", "CalibrationSoftware", "pytsf-demo"),
        ("+", "MzStandardDeviationPPM", "0.5"),
        ("+", "ReferenceMassList", "synthetic-demo"),
        ("+", "ReferencePeakMasses",
         struct.pack(f"<{ref_mz.size}d", *ref_mz)),
        ("+", "MeasuredTimesOfFlight",
         struct.pack(f"<{ref_t.size}d", *ref_t)),
    ])
    con.commit()
    con.close()

    size = sum(os.path.getsize(os.path.join(out, f)) for f in os.listdir(out))
    print(f"wrote {out}: {len(frames)} spectra ({WIDTH}x{HEIGHT}), "
          f"{size / 1e6:.1f} MB")
    print(f"synthetic ions with spatial structure: "
          f"{', '.join(f'{m:.4f}' for m in sorted(patterns))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
