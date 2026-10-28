#!/usr/bin/env python3
"""Generate a small synthetic MALDI-TIMS imaging TDF dataset.

The TDF notebook uses this so that the worked example is self-contained and
reproducible without anybody's real data. The output is a valid `.d` directory
(``analysis.tdf`` + ``analysis.tdf_bin``) that pytsf reads. The documented
``MzCalibration`` / ``TimsCalibration`` columns are present holding placeholder
values; pytsf does not read them.

    python make_demo_tdf_dataset.py [output.d]

What is in it (all invented):

* a 24 x 16 pixel grid, one frame per pixel, 400 mobility scans per frame;
* a tune-mix-like calibrant present in every pixel, recorded in
  ``CalibrationInfo`` the way a real acquisition records it, so that
  ``TdfData``'s automatic mobility calibration has something to find;
* three "sample" ions with different spatial patterns and mobilities, one of
  which is an *isobaric pair*: two species at the same m/z, different 1/K0,
  different spatial distributions -- the case where ion mobility earns its keep;
* a flat matrix background.

It is synthetic: peak positions, mobilities and spatial patterns are made up,
and so is the calibration. It is for exercising the code path, not for science.
"""

from __future__ import annotations

import os
import shutil
import sqlite3
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pytsf.tdf_bin import encode_frame_payload  # noqa: E402

SCALE = 1e12


def published_mz_to_time(ml1, ml2, ml3, m):
    """The published three-point TOF equation (OpenMS TOFCalibration; Gobom et
    al., Anal. Chem. 2002, 74(15), 3915-3923).
    """
    m = np.asarray(m, dtype=float)
    return ml2 + np.sqrt(SCALE / ml1 * m) + ml3 * m


SCHEMA = """
CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
-- pytsf reads only Id / DigitizerTimebase / DigitizerDelay from MzCalibration
-- and nothing from TimsCalibration; the other documented columns hold 1.0.
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
CREATE TABLE MaldiFrameInfo (
    Frame INTEGER PRIMARY KEY, Chip INTEGER, SpotName TEXT,
    RegionNumber INTEGER, XIndexPos INTEGER, YIndexPos INTEGER,
    LaserPower REAL, NumLaserShots INTEGER, LaserRepRate REAL,
    MotorPositionX REAL, MotorPositionY REAL, MotorPositionZ REAL);
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

# A plausible synthetic m/z calibration, in the published equation's own terms.
CAL = {"timebase": 0.125, "delay": 44000.0, "ml1": 152000.0, "ml2": 44000.0 + 18000.0, "ml3": -4.0e-4}

WIDTH, HEIGHT = 24, 16
NUM_SCANS = 400
NUM_SAMPLES = 900000
K0_HI, K0_LO = 1.60, 0.60  # 1/K0 at scan 0 and at the last scan (linear)

# The calibrant: name, m/z, 1/K0 -- values in the spirit of a tune mix, invented.
CALIBRANT = [
    ("cal-322", 322.0481, 0.7363),
    ("cal-622", 622.0290, 0.9915),
    ("cal-922", 922.0098, 1.1986),
    ("cal-1222", 1221.9906, 1.3934),
]

# Sample ions: m/z, 1/K0, spatial pattern name. 760.5851 appears twice with
# different mobility and different distribution: the isobaric pair.
SAMPLE_IONS = [
    (496.3398, 1.030, "gradient"),
    (760.5851, 1.180, "blobs"),
    (760.5851, 1.305, "ring"),
    (885.5499, 1.260, "stripe"),
]


def spatial_patterns(w: int, h: int) -> dict[str, np.ndarray]:
    yy, xx = np.mgrid[0:h, 0:w]
    xn, yn = xx / (w - 1), yy / (h - 1)
    grad = 0.2 + 0.8 * (xn + yn) / 2.0
    blobs = np.exp(-(((xn - 0.3) ** 2 + (yn - 0.35) ** 2) / 0.02)) + np.exp(
        -(((xn - 0.75) ** 2 + (yn - 0.7) ** 2) / 0.012)
    )
    blobs = 0.05 + 0.95 * blobs / blobs.max()
    r = np.sqrt((xn - 0.5) ** 2 + (yn - 0.5) ** 2)
    ring = np.exp(-((r - 0.32) ** 2) / 0.006)
    ring = 0.05 + 0.95 * ring / ring.max()
    stripe = 0.1 + 0.9 * (np.sin(6.0 * xn) > 0.3)
    return {"gradient": grad, "blobs": blobs, "ring": ring, "stripe": stripe}


def main(argv: list[str]) -> int:
    out = argv[1] if len(argv) > 1 else "demo_tims_imaging.d"
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)

    rng = np.random.default_rng(20261020)

    def mz_to_index(m):
        t = published_mz_to_time(CAL["ml1"], CAL["ml2"], CAL["ml3"], m)
        return (t - CAL["delay"]) / CAL["timebase"]

    def k0_to_scan(k0):
        return (k0 - K0_HI) / (K0_LO - K0_HI) * (NUM_SCANS - 1)

    patterns = spatial_patterns(WIDTH, HEIGHT)
    background_mz = np.arange(330.0, 1000.0, 23.7)
    background_k0 = 0.75 + (background_mz - 330.0) / (1000.0 - 330.0) * 0.6

    def add_ion(scan, tof, inten, mz, k0, height, width_scans=4.0):
        """One ion as a mobility peak (Gaussian in scan) at a fixed TOF index,
        spread over the two neighbouring TOF bins so its centroid is exact.
        """
        x = float(mz_to_index(mz))
        x0 = int(np.floor(x))
        frac = x - x0
        s0 = k0_to_scan(k0)
        for s in range(int(s0 - 3 * width_scans), int(s0 + 3 * width_scans) + 1):
            if not 0 <= s < NUM_SCANS:
                continue
            w = np.exp(-0.5 * ((s - s0) / width_scans) ** 2)
            for xi, wi in ((x0, 1.0 - frac), (x0 + 1, frac)):
                v = int(round(height * w * wi))
                if v > 0:
                    scan.append(s)
                    tof.append(xi)
                    inten.append(v)

    blob = bytearray()
    frames, maldi = [], []
    fid = 0
    for y in range(HEIGHT):
        for x in range(WIDTH):
            fid += 1
            scan, tof, inten = [], [], []
            for _, cmz, ck0 in CALIBRANT:
                add_ion(scan, tof, inten, cmz, ck0, 3000.0 * rng.uniform(0.9, 1.1))
            for bmz, bk0 in zip(background_mz, background_k0):
                add_ion(
                    scan,
                    tof,
                    inten,
                    bmz + rng.normal(0, 0.002),
                    bk0,
                    rng.uniform(60, 200),
                    width_scans=5.0,
                )
            for smz, sk0, pat in SAMPLE_IONS:
                h = 150.0 + 4000.0 * float(patterns[pat][y, x]) * rng.uniform(0.85, 1.15)
                add_ion(scan, tof, inten, smz + rng.normal(0, 0.002), sk0, h)
            # a sprinkle of random noise peaks
            n_noise = 300
            scan += list(rng.integers(0, NUM_SCANS, n_noise))
            tof += list(rng.integers(int(mz_to_index(300.0)), int(mz_to_index(1300.0)), n_noise))
            inten += list(rng.integers(1, 25, n_noise))

            scan = np.asarray(scan)
            tof = np.asarray(tof, dtype=np.int64)
            inten = np.asarray(inten, dtype=np.int64)
            # sort by scan then TOF, and merge duplicates (same scan, same TOF)
            order = np.lexsort((tof, scan))
            scan, tof, inten = scan[order], tof[order], inten[order]
            key = scan.astype(np.int64) * (NUM_SAMPLES + 1) + tof
            _uniq, first = np.unique(key, return_index=True)
            summed = np.add.reduceat(inten, first)
            scan, tof, inten = scan[first], tof[first], summed.astype(np.uint32)

            offset = len(blob)
            blob += encode_frame_payload(NUM_SCANS, scan, tof, inten)
            frames.append(
                (
                    fid,
                    fid * 0.1,
                    offset,
                    int(inten.max()),
                    int(inten.sum()),
                    NUM_SCANS,
                    int(inten.size),
                )
            )
            maldi.append(
                (
                    fid,
                    0,
                    f"R00X{x:03d}Y{y:03d}",
                    0,
                    x + 1,
                    y + 1,
                    30.0,
                    200,
                    10000.0,
                    20.0 * x,
                    -20.0 * y,
                    0.0,
                )
            )

    with open(os.path.join(out, "analysis.tdf_bin"), "wb") as f:
        f.write(bytes(blob))

    con = sqlite3.connect(os.path.join(out, "analysis.tdf"))
    con.executescript(SCHEMA)
    con.executemany(
        "INSERT INTO GlobalMetadata VALUES (?,?)",
        {
            "SchemaType": "TDF",
            "SchemaVersionMajor": "3",
            "SchemaVersionMinor": "8",
            "TimsCompressionType": "2",
            "DigitizerNumSamples": str(NUM_SAMPLES),
            "PeakListIndexScaleFactor": "1",
            "MaxNumPeaksPerScan": "500",
            "MzAcqRangeLower": "300.000000",
            "MzAcqRangeUpper": "1300.000000",
            "OneOverK0AcqRangeLower": str(K0_LO),
            "OneOverK0AcqRangeUpper": str(K0_HI),
            "AcquisitionSoftware": "pytsf-demo",
            "AcquisitionSoftwareVersion": "0",
            "InstrumentName": "synthetic",
            "MaldiApplicationType": "Imaging",
        }.items(),
    )
    con.execute(
        "INSERT INTO MzCalibration (Id,ModelType,DigitizerTimebase,DigitizerDelay,"
        "T1,T2,dC1,dC2,C0,C1,C2,C3,C4,C5,C6,C7) VALUES (1,1,?,?,1,1,1,1,1,1,1,1,1,1,1,1)",
        (CAL["timebase"], CAL["delay"]),
    )
    con.execute("INSERT INTO TimsCalibration VALUES (1,1,1,1,1,1,1,1,1,1,1,1)")
    con.executemany(
        "INSERT INTO Frames (Id,Time,Polarity,ScanMode,MsMsType,TimsId,MaxIntensity,"
        "SummedIntensities,NumScans,NumPeaks,MzCalibration,T1,T2,TimsCalibration,"
        "PropertyGroup,AccumulationTime,RampTime,Pressure) "
        "VALUES (?,?,'+',8,0,?,?,?,?,?,1,25.0,23.0,1,NULL,100.0,100.0,2.4)",
        frames,
    )
    con.executemany("INSERT INTO MaldiFrameInfo VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", maldi)

    # The calibrant, recorded as a real acquisition records it: known masses,
    # the flight times they were measured at, and their known mobilities.
    cmz = np.array([c[1] for c in CALIBRANT])
    ck0 = np.array([c[2] for c in CALIBRANT])
    ct = published_mz_to_time(CAL["ml1"], CAL["ml2"], CAL["ml3"], cmz)
    con.executemany(
        "INSERT INTO CalibrationInfo VALUES (?,?,?)",
        [
            ("+", "CalibrationSoftware", "pytsf-demo"),
            ("+", "ReferenceMassList", "synthetic-demo"),
            ("+", "ReferenceMobilityList", "synthetic-demo"),
            ("+", "ReferencePeakMasses", struct.pack(f"<{cmz.size}d", *cmz)),
            ("+", "MeasuredTimesOfFlight", struct.pack(f"<{ct.size}d", *ct)),
            ("+", "ReferencePeakMobilities", struct.pack(f"<{ck0.size}d", *ck0)),
            (
                "+",
                "ReferenceMobilityPeakNames",
                b"".join(c[0].encode() + b"\x00" for c in CALIBRANT),
            ),
        ],
    )
    con.commit()
    con.close()

    print(
        f"wrote {out}: {len(frames)} frames ({WIDTH} x {HEIGHT} pixels), "
        f"{NUM_SCANS} scans each, {os.path.getsize(os.path.join(out, 'analysis.tdf_bin')) / 1e6:.1f} MB"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
