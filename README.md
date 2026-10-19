# pytsf

**Read Bruker timsTOF TSF (`.d`) data on macOS — in pure Python, with no vendor library.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)

<!-- TODO before publishing: replace OWNER/REPO throughout with your GitHub path. -->

Bruker's `tdf-sdk` ships `timsdata.dll` for Windows and `libtimsdata.so` for Linux x86-64,
and nothing else. There is no macOS build — Intel or Apple Silicon — so the standard route
from a timsTOF `.d` folder to anything else simply does not run on a Mac.

`pytsf` reimplements the TSF reading path in Python. It needs NumPy and `zstandard`, and runs
anywhere they do.

```python
from pytsf import TsfData

with TsfData("analysis.d") as tsf:
    print(len(tsf), "spectra")
    mz, intensity = tsf.line_spectrum_mz(1)       # centroided
    mz, intensity = tsf.profile_spectrum_mz(1)    # profile
```

## Install

```bash
pip install git+https://github.com/OWNER/REPO.git
```

or, from a clone:

```bash
pip install -e .
```

## What it does

- **Line (centroided) spectra**, **profile spectra**, and **peak widths**
- **index ↔ m/z** conversion, fitting the published time-of-flight calibration equation to the
  calibration points the instrument recorded in the file (see below)
- The full `analysis.tsf` metadata, plus vendor extension tables such as `MaldiFrameInfo`
  that imaging acquisitions carry
- A drop-in API: method names mirror the Python wrapper in Bruker's own SDK examples, so most
  scripts port by changing the import

**Out of scope:** TDF, the ion-mobility sibling format — use
[OpenTIMS](https://github.com/michalsta/opentims) or
[timsrust](https://github.com/MannLabs/timsrust), which handle it well. Post-hoc recalibrated
states are also unsupported, and raise rather than silently returning raw-calibration values.

## Calling it from C, R, MATLAB, Julia, C#

If you have software that already loads Bruker's shared library by name, `capi/` builds a
replacement exposing the same C functions (`tsf_open`, `tsf_read_line_spectrum_v2`,
`tsf_index_to_mz`, …), so that software can run on macOS. Spectra come back identical; m/z
is computed independently and agrees to about 1 ppm, so it is a drop-in for the *interface*
but not bit-for-bit for the *numbers*:

```bash
cd capi
brew install zstd          # macOS; SQLite ships with the system
make                       # -> libtimsdata_tsf.dylib
make check DATASET=/path/to/analysis.d
```

| target | what it does |
|---|---|
| `make` | build for this machine (the usual case) |
| `make info` | print what the build detects, and whether zstd and SQLite link |
| `make check DATASET=...` | build the smoke test and run it against a dataset |
| `make install PREFIX=...` | install the library and header, default `/usr/local` |
| `make compat-name` | add a `libtimsdata.dylib` symlink for software that loads the vendor library by name |
| `make universal` | macOS arm64 + x86_64 in one binary — needs a universal libzstd, see the notes in the Makefile |

Dependencies are checked before compiling, so a missing one gives an actionable message
rather than a linker error. `make universal` is rarely what you want: Homebrew installs zstd
for your native architecture only, so the universal link fails unless you supply a universal
libzstd yourself. Building separately on each architecture and `lipo`-ing the results is
usually easier.


It is plain C with no Python dependency, needing only SQLite and zstd, and builds on Linux
(`.so`) and macOS (`.dylib`). `capi/smoke_test.c` is a worked example of calling it.

It computes m/z the same way the Python package does, from the same published equation, so a
tool that loads it in place of the vendor library gets ~1 ppm agreement rather than exactness.
That and the other deliberate differences — recalibrated states unsupported,
`tsf_set_num_threads` accepted and ignored, `length <= 0` treated as a pure size query — are
stated at the top of `capi/tsf_compat.h`.

`capi/validate_capi.py` loads this library and the vendor's through the *same* ctypes
bindings and compares every entry point spectrum by spectrum. On both real datasets below,
spectra decode identically through the C ABI and m/z agrees within tolerance inside the
calibrant range, exactly as for the Python package.

## Converting imaging data to imzML

`examples/tsf_to_imzml.ipynb` is a complete notebook: open a MALDI imaging `.d`, convert it to
an imzML/ibd pair with [pyimzML](https://github.com/alexandrovteam/pyimzML), verify the result
round-trips, and render ion images from the converted file. The output opens in Cardinal,
SCiLS, METASPACE, MSiReader and anything else that reads imzML.

It runs out of the box against a small synthetic dataset:

```bash
python examples/make_demo_dataset.py demo_imaging.d
jupyter lab examples/tsf_to_imzml.ipynb
```

Point `ANALYSIS_DIR` at a real `.d` folder to use it for real. A 2,400-pixel acquisition
converts in about a second.

## Accuracy

Correctness is established by **differential testing**: on a Linux machine with the vendor SDK
present, `pytsf` and Bruker's own library are given identical input and their outputs compared.

Validated against two timsTOF fleX MALDI imaging acquisitions — one per polarity, under
different acquisition software versions and digitizer configurations:

| | positive | negative |
|---|---|---|
| spectra | 3,050 | 2,408 |
| line indices + intensities | 21,447,884 — **0 mismatches** | 13,890,740 — **0 mismatches** |
| profile points | 1,215,943,500 — **0 mismatches** | 1,336,139,000 — **0 mismatches** |
| m/z, within the calibrant range | 0.87 ppm median | 0.48 ppm median |

**35 million peaks and 2.55 billion profile points, zero mismatches.** Spectra decode
bit-exactly. m/z is computed independently from the published equation rather than copied from
the vendor, so it agrees to about 1 ppm rather than exactly — see below. Reading and
converting a 3,000-spectrum run takes about 1.4 times as long as the vendor library does
(roughly 15 s against 11 s).

`tests/validate_real_dataset.py` reproduces this check against any `.d` directory: spectra
must decode identically, and m/z must agree within 5 ppm inside the calibrant range, or it
exits non-zero. Worth running once for a new instrument or acquisition mode.

## How m/z is calculated

A TSF file records, in the `CalibrationInfo` table Bruker documents, the known masses of the
calibrant used (`ReferencePeakMasses`) and the flight times at which those peaks were measured
(`MeasuredTimesOfFlight`). `pytsf` fits the **published three-point time-of-flight calibration
equation**

> **t = ml2 + √(10¹² / ml1 · m) + ml3 · m**

to those points by least squares, converting index to flight time as
`t = DigitizerDelay + index × DigitizerTimebase`.

That equation is neither ours nor Bruker's. It is the published form implemented by
[OpenMS](https://github.com/OpenMS/OpenMS)'s
[`TOFCalibration`](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/classOpenMS_1_1TOFCalibration.html),
whose documentation states that "the procedure is very similar to the one described in"

> **Gobom, J., Mueller, M., Egelhofer, V., Theiss, D., Lehrach, H., Nordhoff, E.
> *A calibration method that simplifies and improves accurate determination of peptide
> molecular masses by MALDI-TOF MS.*
> Analytical Chemistry 2002, 74 (15), 3915–3923.**

So the m/z axis rests on exactly three things: that published equation, the columns Bruker
documents in the schema reference distributed with their SDK, and the contents of your own
file. Nothing else.

### How accurate is it?

This is an **external calibration**: accurate across the mass range the calibrant covered,
degrading when extrapolating beyond it, exactly as any external calibration does. Measured
against the vendor library on the two real datasets:

| | positive | negative |
|---|---|---|
| within the calibrant range | **0.87 ppm** median, 2.2 max | **0.48 ppm** median, 0.70 max |
| extrapolated beyond it | 11.1 ppm median, 35 max | 1.2 ppm median, 2.8 max |
| *instrument's own calibration SD* | *0.56 ppm* | *0.091 ppm* |
| [`timsrust-tsf`](https://github.com/MannLabs/timsrust) for comparison | 16.9 ppm median | 1.8 ppm median |

`index_to_mz` warns when asked to extrapolate outside the calibrant range, and returns `NaN`
for flight times before the fitted time origin, where the equation has no physical solution.

For comparison, `timsrust-tsf` interpolates the m/z axis between the acquisition-range
endpoints instead of using the recorded calibration points; its own README states ±10–20 ppm
for its pure-Rust path, which matches what we measure.

## Limitations

Please read these before relying on the output.

- Validation covers **two timsTOF fleX MALDI imaging acquisitions**, both polarities. It does
  **not** cover non-imaging TSF acquisitions, other instrument models, peak widths (absent
  from both files, so that path is tested only synthetically), or recalibrated states.
- The calibration is **one fit per analysis**, and accuracy **degrades outside the
  calibrant's mass range** — quantified above, and `index_to_mz` warns when it happens. If
  your acquisition range extends well below your lowest calibrant mass, as the positive-mode
  dataset's does, the low-m/z end will be the least accurate; the fix is to calibrate over the
  range you acquire. If a file records no calibration points, references more than one
  calibration row, or mixes polarities, `pytsf` refuses at open time rather than guessing.
- The Gauss-Newton fit converges slightly differently in C and in Python, giving m/z values
  that differ by up to ~0.01 ppm between them — a hundredfold below the method's own accuracy,
  but not bit-identical. Spectrum decoding *is* bit-identical in both.
- Intensities are stored as 32-bit floats, so very large counts come back rounded. That is the
  storage format, not a reader limitation.

If your data differs from the validated configurations in any of these respects, run the
differential check against your own files before trusting the output.

## Tests

```bash
pip install -e ".[test]"

pytest tests/                                        # self-consistency only
TIMSDATA_SO=/path/to/libtimsdata.so pytest tests/    # + differential tests

# against a real dataset (the check that matters):
TIMSDATA_SO=/path/to/libtimsdata.so \
    python tests/validate_real_dataset.py /path/to/analysis.d
```

The differential tests skip automatically when the vendor library is absent, which is the
normal situation on macOS. The vendor library is used **only** as a test oracle and is never
imported by the package itself.

## Provenance and acknowledgements

`pytsf` contains **no Bruker source code and redistributes no Bruker binaries**. It was not
produced by decompiling or disassembling the vendor library. Its understanding of the format
comes from:

- **[timsrust](https://github.com/MannLabs/timsrust)** (MannLabs) — in particular the
  `timsrust-tsf` crate, which to our knowledge is the first open-source TSF reader and which
  documents the structure of the binary spectrum file. `pytsf`'s binary reader builds directly
  on it; the applicable licence notice is reproduced in [LICENSE](LICENSE), which also records an inconsistency in upstream's own licence metadata. If you work in Rust,
  use timsrust rather than this.
- **Bruker's own distributed documentation** — the schema reference and public C API header
  shipped with the tdf-sdk, the former of which states that it is the primary reference for
  the file format.
- **Published mass-spectrometry literature** on time-of-flight calibration, in particular the
  calibration model described in [OpenMS](https://github.com/OpenMS/OpenMS)'s
  `TOFCalibration` (after Gobom et al., *Anal. Chem.* 2002).
- **Behavioural testing** through the vendor library's documented public C API.

Related projects worth knowing about:

- **[OpenTIMS / TimsPy](https://github.com/michalsta/opentims)** — the reference open-source
  TDF reader, and the paper that opened up timsTOF data generally
  ([Pietrzak et al., *J. Proteome Res.* 2021](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00962)).
- **[AlphaTims](https://github.com/MannLabs/alphatims)** — fast TDF access for proteomics.
- **[pyimzML](https://github.com/alexandrovteam/pyimzML)** — imzML reading and writing, used
  by the conversion notebook.
- **[OpenMS](https://github.com/OpenMS/OpenMS)** — TSF support is tracked in
  [issue #10552](https://github.com/OpenMS/OpenMS/issues/10552).

This project documents **no detail of the format beyond what the sources above already make
public**. Where behaviour had to be determined empirically, the result is expressed in the
implementation and validated against the vendor library, not published as a specification.

### On the vendor SDK licence

Bruker's tdf-sdk end-user licence prohibits decompiling, disassembling, reverse engineering or
modifying the software, and `pytsf` was written to respect that. The vendor library is
exercised only through the public C API it was published for, and only in the test suite.

Users should be aware that interoperability exceptions to such terms exist in some
jurisdictions — EU Directive 2009/24/EC Art. 6, and more narrowly US DMCA §1201(f) — but this
project does not rely on them. **None of this is legal advice.** If you intend to redistribute
`pytsf` or build a product on it, take your own.

## Citing

If `pytsf` is useful in published work, please cite this repository
(`https://github.com/OWNER/REPO`), and please also cite
[timsrust](https://github.com/MannLabs/timsrust) and the
[OpenTIMS paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00962) — this would not
exist without that groundwork.

## Contributing

Issues and pull requests welcome. The single most valuable contribution is **a TSF dataset
from a configuration not yet covered** — a non-imaging acquisition, a non-fleX instrument, a
file with peak widths, or one whose calibrant does not cover its acquisition range. A small file is
enough.

## Licence

MIT — see [LICENSE](LICENSE), which also carries the third-party notice for
[timsrust](https://github.com/MannLabs/timsrust), from which the binary reader derives.
