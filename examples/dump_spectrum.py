#!/usr/bin/env python3
"""Print a summary of a TSF analysis and dump one spectrum.

This is the pytsf equivalent of the SDK's examples/py/tsfdata_example.py, so you
can compare the two side by side on a machine that has the Bruker library.

    python examples/dump_spectrum.py /path/to/analysis.d [spectrum_id]
"""

from __future__ import annotations

import sys

import numpy as np

from pytsf import TsfData


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2

    analysis_dir = argv[1]
    spectrum_id = int(argv[2]) if len(argv) > 2 else 1

    with TsfData(analysis_dir) as tsf:
        meta = tsf.meta
        print(f"Analysis has {len(tsf)} spectra.")
        print(
            f"Line spectra: {meta.has_line_spectra} "
            f"(peak width: {meta.has_line_spectra_peak_width}); "
            f"profile spectra: {meta.has_profile_spectra}"
        )
        lo, hi = meta.mz_acq_range
        print(f"m/z acquisition range: {lo} .. {hi}")
        print(f"Digitizer samples: {meta.digitizer_num_samples}")
        if meta.is_maldi:
            print("MALDI imaging analysis (MaldiFrameInfo present).")

        print(f"Calibration: {tsf.calibration}")

        frame = meta.frame(spectrum_id)
        print(
            f"\nSpectrum {frame.id}: t={frame.time:.3f}s polarity={frame.polarity} "
            f"MsMsType={frame.msms_type} NumPeaks={frame.num_peaks} "
            f"TIC={frame.summed_intensities} BPI={frame.max_intensity}"
        )

        if not frame.has_data:
            print("  (no data stored for this spectrum)")
            return 0

        if meta.has_line_spectra:
            if meta.has_line_spectra_peak_width:
                idx, inten, width = tsf.read_line_spectrum_with_width(spectrum_id)
                mz = tsf.index_to_mz(spectrum_id, idx)
                print(f"\nLine spectrum: {len(mz)} peaks (with widths)")
                for i in range(min(10, len(mz))):
                    print(f"  mz={mz[i]:14.6f}  I={inten[i]:12.1f}  w={width[i]:.5f}")
            else:
                idx, inten = tsf.read_line_spectrum(spectrum_id)
                mz = tsf.index_to_mz(spectrum_id, idx)
                print(f"\nLine spectrum: {len(mz)} peaks")
                for i in range(min(10, len(mz))):
                    print(f"  mz={mz[i]:14.6f}  I={inten[i]:12.1f}")
            if len(mz):
                print(f"  summed intensity = {float(np.sum(inten)):.1f}")

        if meta.has_profile_spectra:
            inten = tsf.read_profile_spectrum(spectrum_id)
            nz = np.nonzero(inten)[0]
            print(f"\nProfile spectrum: {inten.size} points, {nz.size} non-zero")
            if nz.size:
                mz = tsf.index_to_mz(spectrum_id, nz[:10].astype(np.float64))
                for i in range(min(10, nz.size)):
                    print(f"  idx={nz[i]:8d}  mz={mz[i]:14.6f}  I={inten[nz[i]]:10d}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
