# Provenance

`pytsf` is an independent implementation of a reader for Bruker's TSF (timsTOF
Spectrum Format) data. This document records where its understanding of the
format comes from, and how the implementation was verified.

## Sources

| Component | Source |
|---|---|
| `analysis.tsf` SQLite schema | the schema reference Bruker distributes with the tdf-sdk, which states that it is the primary reference for the file format |
| Public API surface and semantics | the public C API header distributed with the SDK |
| `analysis.tsf_bin` structure | [`timsrust-tsf`](https://github.com/MannLabs/timsrust) (MannLabs; MIT per its LICENSE file — see our LICENSE for a note on upstream's licence metadata) |
| m/z calibration | the three-point TOF equation `t = ml2 + sqrt(1e12/ml1*m) + ml3*m` published in [OpenMS](https://github.com/OpenMS/OpenMS)'s `TOFCalibration`, which cites Gobom *et al.*, *Analytical Chemistry* 2002, **74** (15), 3915–3923; fitted to the reference masses and measured flight times the file records in `CalibrationInfo` — a table documented in Bruker's own schema reference |
| Everything else | behavioural testing through the vendor library's documented public C API |

## What this project does not contain

- **No Bruker source code.** None was available; none was used.
- **No Bruker binaries.** Neither `timsdata.dll` nor `libtimsdata.so` is bundled
  or required at runtime. `pytsf` depends only on NumPy and `zstandard`.
- **No decompiled or disassembled output.** The implementation was not obtained
  by decompiling, disassembling, or inspecting debug symbols of the vendor
  library.
- **No format specification document.** Where behaviour had to be determined
  empirically, the result is expressed in the implementation and validated
  against the vendor library, rather than written up as a specification. Note
  that this is a difference of form, not a claim of secrecy: working code
  necessarily embodies what it implements. The implementation does not
  reproduce the vendor's own calibration algorithm at all; see the note below.

## A note on what this package deliberately does not do

Bruker's schema reference documents that each file carries a `MzCalibration`
row holding a model type and a set of coefficient columns, and that the
spectrum temperatures are "required to perform a temperature compensated mz
calibration". How the vendor library combines those values into an m/z axis
is not published, and Bruker's SDK documentation describes its calibration
algorithms as proprietary. `timsrust` declines to reimplement them for that
reason.

`pytsf` does not reimplement them either. It calculates m/z only from the
published Gobom/OpenMS equation fitted to the calibration points the instrument
itself recorded in the documented `CalibrationInfo` table, and reads from
`MzCalibration` only the two digitiser-timing columns whose names carry their
ordinary meaning. The cost is accuracy: about 1 ppm across the calibrant range
rather than exact agreement with the vendor. That trade is deliberate, and the
README quantifies it.

## Role of the vendor library

Bruker's library is used **only as a test oracle**, and only in `tests/` — never
by the `pytsf` package itself. On a Linux x86-64 machine with a licensed copy of
the SDK, the test suite feeds identical inputs to `pytsf` and to the vendor
library through its documented public C API and asserts that the outputs agree.
This is differential testing against a reference implementation, the ordinary
way to validate a reimplementation, and it establishes correctness without
inspecting the library's internals.

Users without the SDK can still run the full suite: the oracle-dependent tests
skip automatically, and the remainder validate `pytsf` against round-trips and
checked-in fixtures.

## Validation on real data

The implementation was confirmed against two Bruker-written datasets: timsTOF
fleX MALDI imaging acquisitions in positive mode (3050 spectra, acquisition
software 7.0.5) and negative mode (2408 spectra, software 6.1.6, a different
digitizer configuration).

* Line and profile spectra decoded bit-exactly across both files: 35,338,624
  peak indices and intensities and 2,552,082,500 profile points, zero
  mismatches against the vendor library.

Not exercised by either dataset, and so still unvalidated against real data:
non-imaging TSF acquisitions, instruments other than the timsTOF fleX, peak
widths, and recalibrated states.

Both checks are reproducible: `tests/validate_real_dataset.py` compares the
Python implementation against the vendor library on any `.d` directory, and
`capi/validate_capi.py` does the same for the C library through the shared C
ABI. On both datasets above the C library decodes spectra identically to the
vendor and to the Python implementation; its m/z agrees with the Python
implementation to within 0.01 ppm (the two Gauss-Newton fits converge slightly
differently) and with the vendor to about 1 ppm inside the calibrant range.

## Prior art

[`timsrust-tsf`](https://github.com/MannLabs/timsrust) is, to our
knowledge, the first open-source TSF reader, and `pytsf`'s binary reader builds
directly on it. `pytsf` extends it in two respects: it reads profile spectra and
peak widths in addition to line spectra, and it fits the published calibration
equation to the reference points the file records rather than interpolating the
m/z axis between the acquisition-range endpoints. The difference that makes is
quantified in the README.

For TDF — the ion-mobility sibling format — use
[OpenTIMS](https://github.com/michalsta/opentims) or
[timsrust](https://github.com/MannLabs/timsrust); `pytsf` does not read TDF.

## On the vendor SDK licence

Bruker's tdf-sdk end-user licence prohibits decompiling, disassembling, reverse
engineering or modifying the software. `pytsf` was written to respect that: its
understanding of the format comes from Bruker's own distributed documentation,
from open-source prior art, and from published literature, and the vendor
library was exercised only through the public API it was published for.

Interoperability exceptions to such terms exist in some jurisdictions — EU
Directive 2009/24/EC Art. 6 and, more narrowly, US DMCA §1201(f) — but this
project does not rely on them. Nothing here is legal advice.
