#!/usr/bin/env python3
"""Print how a TDF analysis was calibrated and dump one frame as a table.

    python examples/dump_tdf_frame.py /path/to/run.d [frame_id] [--csv out.csv]

Columns: scan, 1/K0, CCS (singly charged), m/z, intensity (as stored).
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
from pytsf import TdfData


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    d = argv[1]
    frame_id = int(argv[2]) if len(argv) > 2 and not argv[2].startswith("--") else None
    csv = argv[argv.index("--csv") + 1] if "--csv" in argv else None

    with warnings.catch_warnings(record=True) as notes:
        warnings.simplefilter("always")
        tdf = TdfData(d)
    for n in notes:
        print(f"note: {n.message}", file=sys.stderr)

    with tdf:
        print(f"{d}: {len(tdf)} frames x {tdf.meta.num_scans} scans")
        print(
            f"m/z   : fitted to {tdf.calibration.ref_mz.size} calibration points, "
            f"{tdf.calibration.mz_min:.1f}-{tdf.calibration.mz_max:.1f} m/z, "
            f"fit RMS {tdf.calibration.fit_rms_ppm:.2f} ppm"
        )
        print(f"1/K0  : {tdf.mobility}")
        for ion in tdf.reference_ions:
            flag = "used" if ion in tdf.mobility.reference_ions else "rejected"
            print(
                f"        {ion.name:>18s}  S/N {ion.snr:8.1f}  m/z err {ion.mz_error_ppm:+6.2f} ppm  "
                f"apex scan {ion.apex_scan:8.2f}  {flag}"
            )

        fid = frame_id if frame_id is not None else tdf.frames[0].id
        mz, k0, inten = tdf.frame_mz_mobility(fid)
        scan = tdf.read_frame(fid).scan
        ccs = tdf.oneoverk0_to_ccs(k0, mz, charge=1)
        print(f"\nframe {fid}: {mz.size} peaks")
        rows = np.column_stack([scan, k0, ccs, mz, inten])
        if csv:
            np.savetxt(
                csv,
                rows,
                delimiter=",",
                header="scan,oneoverk0,ccs_z1,mz,intensity",
                fmt=["%d", "%.5f", "%.2f", "%.5f", "%d"],
                comments="",
            )
            print(f"wrote {csv}")
        else:
            print(f"{'scan':>6s} {'1/K0':>8s} {'CCS(z=1)':>9s} {'m/z':>11s} {'intensity':>10s}")
            for r in rows[:40]:
                print(f"{int(r[0]):6d} {r[1]:8.4f} {r[2]:9.2f} {r[3]:11.5f} {int(r[4]):10d}")
            if rows.shape[0] > 40:
                print(f"... {rows.shape[0] - 40} more; use --csv to get them all")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
