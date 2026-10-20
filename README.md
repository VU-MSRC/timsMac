# pytsf and timsMac

**Read Bruker timsTOF `.d` data — TSF and TDF (ion mobility) — on macOS, in pure Python, with no vendor library.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)

<!-- TODO before publishing: replace OWNER/REPO throughout with your GitHub path. -->

Bruker's `tdf-sdk` ships `timsdata.dll` for Windows and `libtimsdata.so` for Linux x86-64,
and nothing else. There is no macOS build — Intel or Apple Silicon — so the standard route
from a timsTOF `.d` folder to anything else simply does not run on a Mac.

`pytsf` reimplements the reading path in Python for both formats: **TSF** (`analysis.tsf`,
MALDI imaging and other non-mobility acquisitions) and **TDF** (`analysis.tdf`, TIMS
ion-mobility data — PASEF, PRM-PASEF, direct infusion). It needs NumPy and `zstandard`, and
runs anywhere they do. **Not using Python?** A pre-built macOS library with Bruker's own C
function names, for MATLAB, R, Julia, C and anything else, is on the
[Releases page]([https://github.com/OWNER/REPO/releases](https://github.com/VU-MSRC/timsMac/) — see
[Calling it from C, R, MATLAB, Julia, C#](#calling-it-from-c-r-matlab-julia-c).

```python
from pytsf import TsfData, TdfData

with TsfData("imaging.d") as tsf:
    mz, intensity = tsf.line_spectrum_mz(1)       # centroided
    mz, intensity = tsf.profile_spectrum_mz(1)    # profile

with TdfData("pasef.d") as tdf:
    mz, oneoverk0, intensity = tdf.frame_mz_mobility(1)   # every peak of frame 1
    ccs = tdf.oneoverk0_to_ccs(oneoverk0, mz, charge=1)
```

## Install

```bash
pip install git+https://github.com/OWNER/REPO.git
```

or, from a clone:

```bash
pip install -e .                 # the package
pip install -e ".[examples]"     # + matplotlib, pyimzML and JupyterLab for the notebooks
```

## What it does

**TSF**

- **Line (centroided) spectra**, **profile spectra**, and **peak widths**
- **index ↔ m/z** conversion, fitting the published time-of-flight calibration equation to the
  calibration points the instrument recorded in the file (see below)
- The full `analysis.tsf` metadata, plus vendor extension tables such as `MaldiFrameInfo`
  that imaging acquisitions carry
- A drop-in API: method names mirror the Python wrapper in Bruker's own SDK examples, so most
  scripts port by changing the import

**TDF (ion mobility)**

- **Frames** decoded to (scan, TOF index, intensity) arrays — centroided peaks per mobility scan
- **index ↔ m/z**, the same fitted calibration as for TSF
- **scan ↔ 1/K0**, calibrated from the reference ions the file records, when they are present
  in the data (see [How mobility is calculated](#how-mobility-is-calculated))
- **1/K0 ↔ CCS** by the Mason-Schamp equation
- The `PasefFrameMsMsInfo`, `Precursors`, `Frames` and other tables from Bruker's documented
  schema, as rows

**Out of scope:** post-hoc recalibrated states (the reader raises rather than silently
returning raw-calibration values), and vendor-side processing such as the SDK's
`tims_extract_centroided_spectrum_for_frame` — you get the stored peaks and do your own
summing. The C library (`capi/`) covers TSF only.

## Calling it from C, R, MATLAB, Julia, C#

If you have software that already loads Bruker's shared library by name, this repository
provides a replacement exposing the same C functions (`tsf_open`, `tsf_read_line_spectrum_v2`,
`tsf_index_to_mz`, …), so that software can run on macOS. Spectra come back identical; m/z
is computed independently and agrees to about 1 ppm, so it is a drop-in for the *interface*
but not bit-for-bit for the *numbers*.

### Use the pre-built library

The [Releases page](https://github.com/VU-MSRC/timsMac/releases) carries a ready-made
`libtimsdata_tsf.dylib` and its header `tsf_compat.h`. The library is **universal** — one
file for Apple Silicon and Intel — runs on **macOS 12 or later**, and depends on nothing but
macOS itself (zstd is compiled in; SQLite comes with the system). 

```bash
cd ~/Downloads
shasum -a 256 libtimsdata_tsf.dylib              # must match the hash on the release page
xattr -d com.apple.quarantine libtimsdata_tsf.dylib   # once; macOS flags browser downloads
sudo cp libtimsdata_tsf.dylib /usr/local/lib/
sudo cp tsf_compat.h /usr/local/include/
```

Then load it from your software by full path, or by name from `/usr/local/lib`:

```python
import ctypes
lib = ctypes.CDLL("/usr/local/lib/libtimsdata_tsf.dylib")   # Python; declare argtypes as in tsf_compat.h
```

```matlab
loadlibrary('/usr/local/lib/libtimsdata_tsf.dylib', 'tsf_compat.h', 'alias', 'tsf')   % MATLAB
```

```bash
cc myprog.c -I/usr/local/include -L/usr/local/lib -ltimsdata_tsf -Wl,-rpath,/usr/local/lib   # C
```

Bruker's own SDK example `tsfdata.py` works with it once its platform check is given a macOS
branch naming this file. For software that asks for the vendor's file name literally,
`sudo ln -s libtimsdata_tsf.dylib /usr/local/lib/libtimsdata.dylib` — local interoperability
only; it is not the vendor library and should not be passed on under that name.

Every function returns an error code on failure and `tsf_get_last_error_string` explains it
in words. The deliberate differences from the vendor library — recalibrated states
unsupported, `tsf_set_num_threads` accepted and ignored, `length <= 0` treated as a pure size
query, ~1 ppm m/z — are stated at the top of `capi/tsf_compat.h`. `capi/smoke_test.c` is a
complete worked example of calling it from C.

### Build it yourself

Only needed if you want to change the source, build for Linux (`.so`), or prefer not to run a
binary someone else compiled. It is plain C with no Python dependency, needing only SQLite
and zstd:

```bash
cd capi
brew install zstd          # macOS; SQLite ships with the system
make                       # -> libtimsdata_tsf.dylib for this machine
make check DATASET=/path/to/analysis.d
```

Or, with no package manager involved at all — download a
[zstd release](https://github.com/facebook/zstd/releases), unpack it, and:

```bash
make ZSTD_SRC=/path/to/zstd-1.5.7     # builds zstd from source, links it statically
```

The library then depends on nothing but macOS itself, which is the best form to hand to
other software, and it avoids Homebrew — useful on macOS versions Homebrew no longer
builds packages for (anything older than the three most recent releases). The pre-built
release binary was made this way, as a universal build:

```bash
export MACOSX_DEPLOYMENT_TARGET=12.0                 # so the result loads on macOS 12+
make universal ZSTD_SRC=/path/to/zstd-1.5.7          # arm64 + x86_64 in one file
```

| target | what it does |
|---|---|
| `make` | build for this machine (the usual case) |
| `make info` | print what the build detects, and whether zstd and SQLite link |
| `make check DATASET=...` | build the smoke test and run it against a dataset |
| `make install PREFIX=...` | install the library and header, default `/usr/local` |
| `make compat-name` | add a `libtimsdata.dylib` symlink for software that loads the vendor library by name |
| `make ZSTD_SRC=…` | build zstd from its source tree and link it statically: no third-party runtime dependency |
| `make universal` | macOS arm64 + x86_64 in one binary — easiest combined with `ZSTD_SRC=`, which builds zstd for both architectures too |

Dependencies are checked before compiling, so a missing one gives an actionable message
rather than a linker error. `make universal` on its own is rarely what you want: Homebrew
installs zstd for your native architecture only, so the universal link fails. Combine it with
`ZSTD_SRC=`, which passes the architecture flags through to the zstd build, or build on each
architecture separately and `lipo` the results. Run `make universal` as the first build in a
fresh tree — a plain `make` beforehand leaves an arm64-only zstd behind.

`capi/validate_capi.py` loads this library and the vendor's through the *same* ctypes
bindings and compares every entry point spectrum by spectrum. On both real datasets below,
spectra decode identically through the C ABI and m/z agrees within tolerance inside the
calibrant range, exactly as for the Python package.

## Reading ion-mobility (TDF) data

```python
from pytsf import TdfData

with TdfData("run.d") as tdf:
    print(len(tdf), "frames,", tdf.meta.num_scans, "scans each")
    print(tdf.mobility)              # how 1/K0 was calibrated, and from what

    fr = tdf.read_frame(1)           # .scan, .tof_indices, .intensities; .scan_offsets
    mz = tdf.index_to_mz(fr.tof_indices)
    k0 = tdf.scan_to_oneoverk0(fr.scan)

    for row in tdf.meta.pasef_msms_info(frame_id=2):   # PASEF isolation windows
        print(row["ScanNumBegin"], row["ScanNumEnd"], row["IsolationMz"])
```

Intensities are returned **as stored in the file**, which is also what timsrust and OpenTIMS
return. The vendor library returns them scaled by a factor that is not documented anywhere
Bruker distributes (1.0 on one of our files, 2.5 on the other, apparently tied to the
accumulation time); relative intensities within a file are unaffected. If you need to compare
absolute intensities with vendor-derived numbers, be aware of this.

### Worked example: `examples/tdf_analysis.ipynb`

A complete notebook for TIMS data: it opens a `.d`, prints the calibration audit trail (which
reference ions were found, their S/N and m/z error, whether they were used or the reader fell
back), draws a frame as an m/z × mobility map, extracts an ion mobilogram and resolves an
isobaric pair with CCS, builds mobility-integrated and mobility-filtered mass spectra, renders
MALDI ion images with and without a mobility window, and exports per-pixel spectra to imzML
with [pyimzML](https://github.com/alexandrovteam/pyimzML) (mobility-integrated, or one
mobility band). It runs out of the box against a synthetic MALDI-TIMS imaging dataset that
contains a calibrant and an isobaric pair:

```bash
python examples/make_demo_tdf_dataset.py demo_tims_imaging.d
jupyter lab examples/tdf_analysis.ipynb
```

For a quick look at a file without the notebook, `examples/dump_tdf_frame.py run.d` prints the
same calibration report and dumps one frame as scan / 1/K0 / CCS / m/z / intensity.

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

**TDF**, validated against a timsTOF fleX direct-infusion tune-mix run and a timsTOF HT
DDA-PASEF blank:

| | tune mix (fleX, 6.1.6) | DDA-PASEF blank (HT, 7.0.5) |
|---|---|---|
| frames × scans | 293 × 2,858 | 2,725 × 1,054 |
| scan structure + TOF indices | 5,138,519 peaks — **0 mismatches** | 8,657,887 peaks — **0 mismatches** |
| intensities | stored values; vendor returns ×2.5 | stored values; vendor returns ×1 (25 peaks differ by 1 count from float rounding) |
| m/z, within the calibrant range | 0.51 ppm median, 1.5 max | 0.68 ppm |
| 1/K0, from reference ions in the data | **0.25 % median, 0.53 % max** | *no calibrant ions in the file* |
| 1/K0, acquisition-range approximation | 2.1 % median, 2.7 max | 3.5 % median, 4.1 max |
| *instrument's own mobility calibration SD* | *0.26 %* | *0.28 %* |

**13.8 million peaks, zero mismatches** in frame decoding. `tests/validate_real_tdf.py`
reproduces this against any TDF `.d` directory.

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

## How mobility is calculated

A TDF file records, in the documented `CalibrationInfo` table, the calibrant used for the
mobility axis: its ions' known masses and reduced mobilities (`ReferencePeakMasses`,
`ReferencePeakMobilities`). If those ions are present in the data — a direct-infusion
tune-mix run, or any acquisition with calibrant in the sample — `pytsf` locates each one by
m/z, finds its apex scan, and fits **1/K0 as a linear function of scan number** through them.
This is the ordinary external calibration of ion-mobility mass spectrometry, and its physical
basis is that TIMS elutes ions by ramping the trapping field linearly in time, with 1/K0 linear
in that field to first order:

> Michelmann, K., Silveira, J. A., Ridgeway, M. E., Park, M. A. *Fundamentals of Trapped Ion
> Mobility Spectrometry.* J. Am. Soc. Mass Spectrom. 2015, 26, 14–24.
>
> Ridgeway, M. E., Lubeck, M., Jordens, J., Mann, M., Park, M. A. *Trapped ion mobility
> spectrometry: A short review.* Int. J. Mass Spectrom. 2018, 425, 22–35.

Reference 1/K0 values for the usual tune mix are tabulated in the CCS compendium of
Stow et al., *Anal. Chem.* 2017, 89, 9048–9055; `pytsf` uses the values the file itself records.

Each located ion must pass quality gates before it is used — signal-to-background ≥ 50,
observed m/z within 2.5 ppm of the known mass, and the final fit residual ≤ 1 % with 1/K0
falling with scan number, as TIMS requires. On a blank these gates reject the noise picks and
the reader **falls back, with a warning,** to the approximation
[timsrust](https://github.com/MannLabs/timsrust) uses: linear interpolation between the
documented `OneOverK0AcqRangeLower/Upper` endpoints. The real voltage ramp overshoots the
nominal range, so that approximation is only good to a few percent (measured above). You can
also pass your own `MobilityCalibration` built from ions you know to be in your data.

CCS is computed from 1/K0 with the Mason-Schamp equation, using CODATA 2018 constants, N₂ as
the buffer gas and the 305 K convention of the TIMS literature; it reproduces the vendor's
CCS values to every displayed digit.

The vendor records its own mobility calibration in a `TimsCalibration` table whose coefficient
layout is not published. `pytsf` does not read it.

## Limitations

Please read these before relying on the output.

- TSF validation covers **two timsTOF fleX MALDI imaging acquisitions**, both polarities. It
  does **not** cover non-imaging TSF acquisitions, other instrument models, peak widths
  (absent from both files, so that path is tested only synthetically), or recalibrated states.
- TDF validation covers **one direct-infusion tune-mix run and one DDA-PASEF blank**. Frame
  decoding is bit-exact on both; it does **not** yet cover MALDI PASEF imaging, DIA-PASEF, or
  files with more than one `MzCalibration` row (refused at open time).
- **Mobility accuracy depends on calibrant ions being in the data.** With them, 1/K0 agrees
  with the vendor to ~0.25 % inside the span those ions cover, degrading when extrapolating
  beyond it (`scan_to_oneoverk0` warns). Without them — a blank, or a sample run without
  calibrant — the acquisition-range approximation is used and 1/K0 is only good to **2–4 %**,
  with a warning. For the best result on real samples, acquire with calibrant present, or pass
  `MobilityCalibration.from_reference_ions(...)` built from ions of known mobility in your
  sample.
- TDF intensities are the stored values, which differ from the vendor library's by an
  undocumented constant factor per file (see above).
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
    python tests/validate_real_dataset.py /path/to/imaging.d     # TSF
TIMSDATA_SO=/path/to/libtimsdata.so \
    python tests/validate_real_tdf.py /path/to/pasef.d           # TDF
```

The differential tests skip automatically when the vendor library is absent, which is the
normal situation on macOS. The vendor library is used **only** as a test oracle and is never
imported by the package itself.

## Provenance and acknowledgements

`pytsf` contains **no Bruker source code and redistributes no Bruker binaries**. It was not
produced by decompiling or disassembling the vendor library. Its understanding of the format
comes from:

- **[timsrust](https://github.com/MannLabs/timsrust)** (MannLabs) — the `timsrust-tsf` crate,
  which to our knowledge is the first open-source TSF reader and documents the structure of
  the binary spectrum file, and the `timsrust-tdf` crate's frame reader, which documents the
  TDF frame layout (as does [OpenTIMS](https://github.com/michalsta/opentims) independently).
  `pytsf`'s binary readers build directly on them; the applicable licence notice is reproduced
  in [LICENSE](LICENSE), which also records an inconsistency in upstream's own licence
  metadata. If you work in Rust, use timsrust rather than this.
- **Bruker's own distributed documentation** — the schema reference and public C API header
  shipped with the tdf-sdk, the former of which states that it is the primary reference for
  the file format.
- **Published mass-spectrometry literature**: the time-of-flight calibration model described
  in [OpenMS](https://github.com/OpenMS/OpenMS)'s `TOFCalibration` (after Gobom et al.,
  *Anal. Chem.* 2002); the TIMS operating principle (Michelmann et al. 2015, Ridgeway et al.
  2018); the Mason-Schamp equation with CODATA constants; and the tune-mix CCS compendium of
  Stow et al. 2017.
- **Behavioural testing** through the vendor library's documented public C API.

Related projects worth knowing about:

- **[OpenTIMS / TimsPy](https://github.com/michalsta/opentims)** — the reference open-source
  TDF reader, and the paper that opened up timsTOF data generally
  ([Pietrzak et al., *J. Proteome Res.* 2021](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00962)).
  For heavy proteomics TDF work it is the more mature choice; `pytsf` adds a macOS-friendly
  pure-Python path and a reference-ion mobility calibration.
- **[AlphaTims](https://github.com/MannLabs/alphatims)** — fast TDF access for proteomics.
- **[pyimzML](https://github.com/alexandrovteam/pyimzML)** — imzML reading and writing, used
  by the conversion notebook.
- **[OpenMS](https://github.com/OpenMS/OpenMS)** — TSF support is tracked in
  [issue #10552](https://github.com/OpenMS/OpenMS/issues/10552).

This project documents **no detail of the format beyond what the sources above already make
public**. Where behaviour had to be determined empirically, the result is expressed in the
implementation and validated against the vendor library.


## Citing

If `pytsf` is useful in published work, please cite this repository
(`https://github.com/OWNER/REPO](https://github.com/VU-MSRC/timsMac/`), and please also cite
[timsrust](https://github.com/MannLabs/timsrust) and the
[OpenTIMS paper](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00962) — this would not
exist without that groundwork.

## Contributing

Issues and pull requests welcome. 

## Licence

MIT — see [LICENSE](LICENSE), which also carries the third-party notice for
[timsrust](https://github.com/MannLabs/timsrust), from which the binary readers derive.
