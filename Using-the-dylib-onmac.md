# Using the pytsf dylib on a Mac


`libtimsdata_tsf.dylib` lets software on any Mac (Intel or Apple Silicon, macOS 12 or later) read Bruker timsTOF **TSF** `.d` data through the same C functions as Bruker's own `timsdata` library, which has no macOS build. This page covers getting the file onto your Mac, where to put it, and calling it from Python, MATLAB, R, Julia and C.

## 1. Get the file

You need two files from the release page, https://github.com/VU-MSRC/timsMac/releases :

| File | What it is |
| --- | --- |
| `libtimsdata_tsf.dylib` | the library: one file, arm64 + x86\_64, macOS 12 or later |
| `tsf_compat.h` | the C header that declares its functions; needed for C, C++ and MATLAB, useful as reference for everything else |

Open Terminal (Applications → Utilities) and, in the folder where the download landed:

```bash
cd ~/Downloads
shasum -a 256 libtimsdata_tsf.dylib
xattr -d com.apple.quarantine libtimsdata_tsf.dylib
file libtimsdata_tsf.dylib
```

- `shasum` prints a 64-character hash; it must equal the one on the release page. For the build of 20 October 2026 it is `db250c621006c889dabe13d3f8692f6da5163982622c4baa754ffd0259669c4f`.
- `xattr -d` removes the quarantine flag macOS puts on anything a browser downloads; without this, a GUI application may refuse to load the library ("cannot be opened because the developer cannot be verified"). If it says `No such xattr`, the flag was not there — fine.
- `file` should say *Mach-O universal binary with 2 architectures* — the same file serves Intel and Apple Silicon Macs.

The library is ad-hoc signed, which is all macOS requires for code a program loads from disk. It is not notarised, so the quarantine step above is needed once per downloaded copy.

## 2. Where to put it

The simplest rule: **load it by its full path.** Every example below accepts a full path, and a full path never depends on search settings.

If you would rather load it by name, put it where macOS looks for libraries:

```bash
sudo mkdir -p /usr/local/lib /usr/local/include
sudo cp libtimsdata_tsf.dylib /usr/local/lib/
sudo cp tsf_compat.h /usr/local/include/
```

macOS searches `/usr/local/lib` by default when a program asks for a bare library name. For a per-user install without `sudo`, use `~/lib` instead, or point the search path at any folder for the current Terminal session:

```bash
export DYLD_LIBRARY_PATH="$HOME/my-libs:$DYLD_LIBRARY_PATH"
```

That variable is not inherited by applications launched from the Dock or Finder — another reason to prefer the full path.

The library's internal install name is `@rpath/libtimsdata_tsf.dylib`, which matters only when linking a C or C++ program against it (section 7): add `-Wl,-rpath,<folder>` so the finished program knows where to look.

## 3. What it provides

Eleven C functions with the same names and signatures as the TSF part of Bruker's `timsdata` API, so code written for that API ports by changing the library it loads. A *spectrum id* is the `Id` column of the `Frames` table in `analysis.tsf`, starting at 1.

| Function | Purpose | Returns |
| --- | --- | --- |
| `tsf_open(dir, 0)` | open a `.d` folder (the folder, not a file inside it); second argument must be 0 | handle, or 0 on failure |
| `tsf_close(handle)` | release it | — |
| `tsf_get_last_error_string(buf, len)` | the last error on this thread, as text | length needed |
| `tsf_read_line_spectrum_v2(handle, id, index[], intensity[], length)` | centroided peaks: index positions (double) and intensities (float) | peak count, or −1 on error |
| `tsf_read_line_spectrum_with_width_v2(…, width[], length)` | the same plus peak widths; fails if the file stores none | peak count, or −1 |
| `tsf_read_profile_spectrum_v2(handle, id, profile[], length)` | the dense profile; array position is the index | point count, or −1 |
| `tsf_index_to_mz(handle, id, in[], out[], count)` | index → m/z | 1 ok, 0 error |
| `tsf_mz_to_index(handle, id, in[], out[], count)` | m/z → index | 1 ok, 0 error |
| `tsf_has_recalibrated_state(handle)` | always 0 (see section 9) | 0 |
| `tsf_set_num_threads(n)` | accepted and ignored | — |
| `tsf_compat_version()` | version string of this library | text |

Three conventions every language section relies on:

1. **Size query first.** Call a read function with `length = 0` (the array pointers may be null): nothing is written and the return value is the number of entries you need. Allocate that many, call again with the real `length`. Peak counts vary per spectrum, so this is the normal pattern, not an edge case.
2. **Index, then m/z.** Spectra come back as *index* positions on the digitiser axis. Pass them through `tsf_index_to_mz` to get m/z. The same calibration covers the whole file, so the spectrum id given to the conversion functions does not change the result.
3. **Errors are a return code plus a message.** A handle of 0, a count of −1, or a conversion result of 0 means failure; `tsf_get_last_error_string` says why in plain words. Call it with `(NULL, 0)` to learn the length, then with a buffer of that size.

## 4. Python

Python's built-in `ctypes` loads the library directly; the only extra package is NumPy. Save this as `read_tsf.py` and run `python3 read_tsf.py /path/to/analysis.d`. (It was run as written against real imaging data; only the library path changes between platforms.)

```python
"""Read a TSF .d folder through libtimsdata_tsf.dylib with ctypes."""
import ctypes, sys
import numpy as np

LIB = "/usr/local/lib/libtimsdata_tsf.dylib"      # full path to the library

lib = ctypes.CDLL(LIB)

# --- declare the signatures (from tsf_compat.h) ------------------------------
lib.tsf_open.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
lib.tsf_open.restype = ctypes.c_uint64
lib.tsf_close.argtypes = [ctypes.c_uint64]
lib.tsf_close.restype = None
lib.tsf_get_last_error_string.argtypes = [ctypes.c_char_p, ctypes.c_uint32]
lib.tsf_get_last_error_string.restype = ctypes.c_uint32
lib.tsf_read_line_spectrum_v2.argtypes = [
    ctypes.c_uint64, ctypes.c_int64,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_float), ctypes.c_int32]
lib.tsf_read_line_spectrum_v2.restype = ctypes.c_int32
lib.tsf_read_profile_spectrum_v2.argtypes = [
    ctypes.c_uint64, ctypes.c_int64, ctypes.POINTER(ctypes.c_uint32), ctypes.c_int32]
lib.tsf_read_profile_spectrum_v2.restype = ctypes.c_int32
conv = [ctypes.c_uint64, ctypes.c_int64,
        ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_uint32]
lib.tsf_index_to_mz.argtypes = conv
lib.tsf_index_to_mz.restype = ctypes.c_uint32
lib.tsf_mz_to_index.argtypes = conv
lib.tsf_mz_to_index.restype = ctypes.c_uint32


def last_error():
    n = lib.tsf_get_last_error_string(None, 0)
    buf = ctypes.create_string_buffer(n)
    lib.tsf_get_last_error_string(buf, n)
    return buf.value.decode(errors="replace")


def read_line_spectrum(handle, spectrum_id):
    """Return (index, intensity) arrays of one centroided spectrum."""
    n = lib.tsf_read_line_spectrum_v2(handle, spectrum_id, None, None, 0)   # size query
    if n < 0:
        raise RuntimeError(last_error())
    index = np.empty(n, dtype=np.float64)
    intensity = np.empty(n, dtype=np.float32)
    if n:
        lib.tsf_read_line_spectrum_v2(
            handle, spectrum_id,
            index.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            intensity.ctypes.data_as(ctypes.POINTER(ctypes.c_float)), n)
    return index, intensity


def read_profile_spectrum(handle, spectrum_id):
    """Return the dense profile (position = index) of one spectrum, or None."""
    n = lib.tsf_read_profile_spectrum_v2(handle, spectrum_id, None, 0)
    if n < 0:
        return None            # e.g. the file holds no profile spectra
    prof = np.empty(n, dtype=np.uint32)
    if n:
        lib.tsf_read_profile_spectrum_v2(
            handle, spectrum_id, prof.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32)), n)
    return prof


def index_to_mz(handle, spectrum_id, index):
    index = np.ascontiguousarray(index, dtype=np.float64)
    mz = np.empty_like(index)
    ok = lib.tsf_index_to_mz(handle, spectrum_id,
                             index.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                             mz.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), index.size)
    if not ok:
        raise RuntimeError(last_error())
    return mz


# --- use it -------------------------------------------------------------------
handle = lib.tsf_open(sys.argv[1].encode("utf-8"), 0)     # 0 = raw calibration
if handle == 0:
    sys.exit("tsf_open failed: " + last_error())

try:
    spectrum_id = 1                                        # Frames.Id, 1-based
    index, intensity = read_line_spectrum(handle, spectrum_id)
    mz = index_to_mz(handle, spectrum_id, index)
    print(f"spectrum {spectrum_id}: {mz.size} peaks")
    for m, i in list(zip(mz, intensity))[:5]:
        print(f"  m/z {m:12.5f}   intensity {i:10.1f}")
    prof = read_profile_spectrum(handle, spectrum_id)
    if prof is not None:
        print(f"  profile: {prof.size} points, max {prof.max()}")
finally:
    lib.tsf_close(handle)
```

Expected output has the form:

```
spectrum 1: 2861 peaks
  m/z    305.52508   intensity       35.0
  m/z    305.72032   intensity       45.0
  ...
  profile: 554875 points, max 10914
```

To loop over a whole acquisition, read the spectrum ids and their pixel positions straight from `analysis.tsf` with Python's built-in `sqlite3` module (`SELECT Id FROM Frames`; for imaging, join `MaldiFrameInfo` for `XIndexPos`, `YIndexPos`).

If you already use Bruker's example `tsfdata.py` from the SDK (`examples/py/`), it needs one edit to run on a Mac: its platform check near the top raises *Unsupported platform* on anything but Windows or Linux. Add a branch for macOS that names this library, for example `elif sys.platform == "darwin": libname = "/usr/local/lib/libtimsdata_tsf.dylib"`. Every `tsf_*` function that file declares is provided here, so the rest of it, and `tsfdata_example.py`, run unchanged.

**Python users have a second option.** The `pytsf` package in the same repository does everything this library does in pure Python, with no `.dylib` at all (`pip install git+https://github.com/OWNER/REPO.git`, then `from pytsf import TsfData`), and additionally reads TDF ion-mobility data. The library exists for software that is not Python, or that already speaks the C API.

## 5. MATLAB

MATLAB's `loadlibrary` reads the C header and exposes every function through `calllib`. It needs a C compiler configured once (`mex -setup C`; Apple's command-line tools, `xcode-select --install`, are enough). *This section follows MATLAB's documented interface and the header; it has not been run under MATLAB.*

MATLAB's header parser does not understand the compiler attribute the header uses for symbol visibility, so make a copy without it first (once, in Terminal):

```bash
sed 's/__attribute__((visibility("default")))//' tsf_compat.h > tsf_compat_matlab.h
```

Then, in MATLAB:

```matlab
lib = '/usr/local/lib/libtimsdata_tsf.dylib';
loadlibrary(lib, 'tsf_compat_matlab.h', 'alias', 'tsf');

h = calllib('tsf', 'tsf_open', '/path/to/analysis.d', uint32(0));
if h == 0
    error('tsf_open failed: %s', tsfLastError());
end

id = int64(1);                                   % Frames.Id
% size query: null pointers, length 0 -> number of peaks
n = calllib('tsf', 'tsf_read_line_spectrum_v2', h, id, ...
            libpointer('doublePtr'), libpointer('singlePtr'), int32(0));
if n < 0, error(tsfLastError()); end

idxP = libpointer('doublePtr', zeros(1, n));
intP = libpointer('singlePtr', zeros(1, n, 'single'));
calllib('tsf', 'tsf_read_line_spectrum_v2', h, id, idxP, intP, int32(n));

mzP = libpointer('doublePtr', zeros(1, n));
calllib('tsf', 'tsf_index_to_mz', h, id, idxP, mzP, uint32(n));

mz = mzP.Value;  intensity = double(intP.Value);
fprintf('spectrum %d: %d peaks, m/z %.4f .. %.4f\n', id, n, min(mz), max(mz));

calllib('tsf', 'tsf_close', h);
unloadlibrary tsf

function msg = tsfLastError()
    len = calllib('tsf', 'tsf_get_last_error_string', libpointer('cstring'), uint32(0));
    buf = libpointer('cstring', blanks(double(len)));
    calllib('tsf', 'tsf_get_last_error_string', buf, uint32(len));
    msg = buf.Value;
end
```

Profile spectra follow the same two-call pattern with `tsf_read_profile_spectrum_v2` and a `uint32Ptr`. An uninitialised `libpointer('doublePtr')` is a null pointer, which is what the size query expects.

If `loadlibrary` complains about the header, the quickest diagnosis is `loadlibrary(lib, 'tsf_compat_matlab.h', 'mfilename', 'tsf_proto')`, which writes the prototypes it understood to `tsf_proto.m` for inspection.

## 6. Julia and R

*Both snippets follow the header's signatures and each language's documented foreign-function interface; neither has been run.*

**Julia** calls C directly with `ccall`; no package needed:

```julia
const LIB = "/usr/local/lib/libtimsdata_tsf.dylib"

function tsf_error()
    n = ccall((:tsf_get_last_error_string, LIB), UInt32, (Ptr{UInt8}, UInt32), C_NULL, 0)
    buf = Vector{UInt8}(undef, n)
    ccall((:tsf_get_last_error_string, LIB), UInt32, (Ptr{UInt8}, UInt32), buf, n)
    return String(buf[1:end-1])
end

h = ccall((:tsf_open, LIB), UInt64, (Cstring, UInt32), "/path/to/analysis.d", 0)
h == 0 && error("tsf_open failed: ", tsf_error())

id = Int64(1)                                             # Frames.Id
n = ccall((:tsf_read_line_spectrum_v2, LIB), Int32,
          (UInt64, Int64, Ptr{Float64}, Ptr{Float32}, Int32), h, id, C_NULL, C_NULL, 0)
n < 0 && error(tsf_error())

index = Vector{Float64}(undef, n); intensity = Vector{Float32}(undef, n)
ccall((:tsf_read_line_spectrum_v2, LIB), Int32,
      (UInt64, Int64, Ptr{Float64}, Ptr{Float32}, Int32), h, id, index, intensity, n)

mz = similar(index)
ccall((:tsf_index_to_mz, LIB), UInt32,
      (UInt64, Int64, Ptr{Float64}, Ptr{Float64}, UInt32), h, id, index, mz, n)

println("spectrum $id: $n peaks, m/z $(minimum(mz)) .. $(maximum(mz))")
ccall((:tsf_close, LIB), Cvoid, (UInt64,), h)
```

**R** has no direct equivalent of `ctypes`; the clean route is a few lines of C++ through Rcpp, compiled against the library once. Save this as `tsf.cpp`:

```cpp
#include <Rcpp.h>
#include <tsf_compat.h>
using namespace Rcpp;

static std::string lastError() {
    uint32_t n = tsf_get_last_error_string(NULL, 0);
    std::string s(n, '\0');
    tsf_get_last_error_string(&s[0], n);
    return s.c_str();
}

// [[Rcpp::export]]
List tsf_read_spectrum(std::string analysis_dir, double spectrum_id) {
    uint64_t h = tsf_open(analysis_dir.c_str(), 0);
    if (!h) stop("tsf_open failed: " + lastError());
    int64_t id = (int64_t)spectrum_id;
    int32_t n = tsf_read_line_spectrum_v2(h, id, NULL, NULL, 0);
    if (n < 0) { std::string e = lastError(); tsf_close(h); stop(e); }
    std::vector<double> index(n); std::vector<float> inten(n);
    tsf_read_line_spectrum_v2(h, id, index.data(), inten.data(), n);
    NumericVector mz(n), intensity(n);
    tsf_index_to_mz(h, id, index.data(), mz.begin(), (uint32_t)n);
    for (int i = 0; i < n; i++) intensity[i] = inten[i];
    tsf_close(h);
    return List::create(_["mz"] = mz, _["intensity"] = intensity);
}
```

Then, in R (`install.packages("Rcpp")` once; the compiler is Apple's command-line tools):

```r
Sys.setenv(PKG_CPPFLAGS = "-I/usr/local/include",
           PKG_LIBS = "-L/usr/local/lib -ltimsdata_tsf -Wl,-rpath,/usr/local/lib")
Rcpp::sourceCpp("tsf.cpp")

s <- tsf_read_spectrum("/path/to/analysis.d", 1)
plot(s$mz, s$intensity, type = "h", xlab = "m/z", ylab = "intensity")
```

For many spectra, keep the handle open across calls rather than reopening per spectrum: hold it in a static variable or an `XPtr`, and expose `tsf_open`/`tsf_close` as their own exported functions. Alternatively, R users with Python available can call the `pytsf` package through `reticulate` and skip the library altogether.

## 7. C and C++

Include the header, link with `-ltimsdata_tsf`, and tell the linker where the library lives at run time with `-rpath`. This program was compiled and run as written (against the same source built as a `.so`):

```c
#include <stdio.h>
#include <stdlib.h>
#include "tsf_compat.h"

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s /path/to/analysis.d\n", argv[0]); return 2; }

    uint64_t h = tsf_open(argv[1], 0);
    if (!h) {
        char msg[512]; tsf_get_last_error_string(msg, sizeof msg);
        fprintf(stderr, "tsf_open: %s\n", msg); return 1;
    }

    int64_t id = 1;                                   /* Frames.Id */
    int32_t n = tsf_read_line_spectrum_v2(h, id, NULL, NULL, 0);   /* size query */
    if (n < 0) { char msg[512]; tsf_get_last_error_string(msg, sizeof msg); fprintf(stderr, "%s\n", msg); return 1; }

    double *index = malloc(n * sizeof *index);
    float  *inten = malloc(n * sizeof *inten);
    double *mz    = malloc(n * sizeof *mz);
    tsf_read_line_spectrum_v2(h, id, index, inten, n);
    tsf_index_to_mz(h, id, index, mz, (uint32_t)n);

    printf("spectrum %lld: %d peaks\n", (long long)id, n);
    for (int i = 0; i < n && i < 5; i++)
        printf("  m/z %12.5f   intensity %10.1f\n", mz[i], inten[i]);

    free(index); free(inten); free(mz);
    tsf_close(h);
    return 0;
}
```

With the library and header installed under `/usr/local` (section 2):

```bash
cc -std=c99 -Wall -o first_peaks first_peaks.c \
   -I/usr/local/include -L/usr/local/lib -ltimsdata_tsf -Wl,-rpath,/usr/local/lib
./first_peaks /path/to/analysis.d
```

With the two files simply sitting next to your source instead, use `-I. -L. -Wl,-rpath,@executable_path` and keep the `.dylib` beside the built program. C++ needs no changes: the header wraps its declarations in `extern "C"`.

Building for both architectures yourself is one flag away, since the library already carries both: add `-arch arm64 -arch x86_64` to the `cc` line. Set `MACOSX_DEPLOYMENT_TARGET=12.0` in the environment first if the program should run on macOS 12-26 as well as 27.

Handles are safe to use from one thread at a time; for parallel reading, open the same `.d` folder once per thread. Error messages are kept per thread.

## 8. Software that loads Bruker's library by name

Some tools do not take a library path; they ask the system for `libtimsdata.dylib` (or `timsdata`) and expect it to be found. For those, give this library that name too, as a symlink, and make sure it is on the search path:

```bash
cd /usr/local/lib
sudo ln -sf libtimsdata_tsf.dylib libtimsdata.dylib
```

Two things to keep in mind when doing that:

- **It is a stand-in for the interface, not the vendor.** Spectra decode identically, but m/z is computed from the published calibration equation and agrees with Bruker's library to about 1 ppm, not exactly. A tool that compares its output against vendor-produced numbers will see that difference. Use the symlink for local interoperability and do not pass the file on under the vendor's name.
- **Only the `tsf_*` functions exist.** A tool that also calls the `tims_*` (TDF) functions will fail to load it, with an error naming the missing symbol. The Python package in the repository reads TDF; this library does not.

Tools that load the library from a fixed folder of their own (next to their executable, or a `libs/` folder as Bruker's Python examples do) take the file there instead; copying it in under whichever name they look for works the same way.

## 9. Things to know

- **TSF only.** `analysis.tsf` + `analysis.tsf_bin`: MALDI imaging and other single-spectrum timsTOF acquisitions. A TDF folder (`analysis.tdf`, ion mobility) is refused at `tsf_open` with a message saying so; use the Python package for those.
- **m/z is an external calibration.** The library fits the published time-of-flight equation (OpenMS `TOFCalibration`, after Gobom et al., *Anal. Chem.* 2002) to the reference masses and flight times the instrument recorded in the file. Inside the calibrant's mass range it agrees with the vendor library to about 1 ppm (measured: 0.5-0.9 ppm median on two imaging files); outside that range it extrapolates and accuracy degrades, as any external calibration does. Flight times before the fitted time origin return `NaN`. The repository README quantifies all of this.
- **One calibration per file.** The spectrum id passed to `tsf_index_to_mz` / `tsf_mz_to_index` is accepted for API compatibility and does not change the result. Files whose spectra reference more than one calibration row are refused at open time rather than guessed at.
- **Recalibrated states are not supported.** `tsf_open(dir, 1)` fails, and `tsf_has_recalibrated_state` always returns 0. The raw acquisition-time calibration is what you get.
- **Intensities are 32-bit floats**, as stored in the file. Very large counts come back rounded; that is the file format.
- **Spectra decode bit-exactly** against the vendor library: validated on 35 million peaks and 2.5 billion profile points across two instruments' files, zero mismatches.

### Errors you may see

| Message | Meaning | What to do |
| --- | --- | --- |
| `cannot be opened because the developer cannot be verified` / `code signature invalid` when loading | quarantine flag from a download | `xattr -d com.apple.quarantine libtimsdata_tsf.dylib` (section 1) |
| `image not found` / `Library not loaded` | the loader cannot find the file | load by full path, or install to `/usr/local/lib` (section 2) |
| `symbol not found: _tims_open` (or any `tims_*`) | the calling software wants the TDF functions | not provided; see section 8 |
| `Cannot open …/analysis.tsf: unable to open database file` | the path is not a TSF `.d` folder, or is a TDF folder | check the path points at the `.d` folder itself and that it contains `analysis.tsf` |
| `Not a TSF analysis: SchemaType=TDF` | ion-mobility data | use the Python package |
| `records no ReferencePeakMasses / MeasuredTimesOfFlight` | the file carries no calibration points, so m/z cannot be computed | nothing in this library can help; the vendor's calibration model is not reimplemented |
| `use_recalibrated_state is not supported` | second argument of `tsf_open` was 1 | pass 0 |
| `contains no profile spectra` / `no peak widths` | the file was acquired without them | read line spectra only |
| `Spectra reference N different MzCalibration rows` | multiple calibrations in one file | not supported by this library |

Every one of these comes back through `tsf_get_last_error_string`, so a program that prints that string on failure shows the reason directly.

Sources: `tsf_compat.h` and `README.md` in the repository (`https://github.com/OWNER/REPO`); the Python and C examples above were executed against real timsTOF fleX imaging data.
