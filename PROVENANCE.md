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


## Prior references

[`timsrust-tsf`](https://github.com/MannLabs/timsrust) is, to our
knowledge, the first open-source TSF reader, and `pytsf`'s binary reader builds
directly on it. `pytsf` extends it in two respects: it reads profile spectra and
peak widths in addition to line spectra, and it fits the published calibration
equation to the reference points the file records rather than interpolating the
m/z axis between the acquisition-range endpoints. The difference that makes is
quantified in the README.


## On the vendor SDK licence

Bruker's tdf-sdk end-user licence prohibits decompiling, disassembling, reverse
engineering or modifying the software. `pytsf` was written to respect that: its
understanding of the format comes from Bruker's own distributed documentation,
from open-source prior art, and from published literature, and the vendor
library was exercised only through the public API it was published for.

Interoperability exceptions to such terms exist in some jurisdictions — EU
Directive 2009/24/EC Art. 6 and, more narrowly, US DMCA §1201(f) — but this
project does not rely on them. Nothing here is legal advice.
