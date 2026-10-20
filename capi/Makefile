# Build the pytsf TSF C-API compatibility library.
#
#   make                 build for this machine          -> libtimsdata_tsf.dylib / .so
#   make universal       macOS: arm64 + x86_64 in one    (needs universal deps, see below)
#   make check           build and run the smoke test against DATASET=...
#   make install         install to $(PREFIX), default /usr/local
#   make compat-name     also create a libtimsdata.dylib symlink (see below)
#   make info            print what the build has detected
#   make clean
#
# ---------------------------------------------------------------------------
# Dependencies
#
#   SQLite  macOS ships it; Linux: apt install libsqlite3-dev
#   zstd    macOS: brew install zstd       Linux: apt install libzstd-dev
#
#           -- or, with no package manager at all --
#
#           make ZSTD_SRC=/path/to/zstd-1.5.x
#
#           builds zstd from its source tree (https://github.com/facebook/zstd,
#           download a release tarball and unpack it) and links it STATICALLY.
#           The resulting library then depends on nothing but the OS itself,
#           which is the best form to hand to other people's software, and it
#           sidesteps Homebrew entirely -- useful on macOS versions Homebrew
#           no longer builds packages for.
#
# ---------------------------------------------------------------------------
# A note on universal builds
#
# `make` produces a library for the architecture you are on, which is what you
# want in almost all cases. `make universal` produces one carrying both arm64
# and x86_64 -- but it can only work if *every* library it links also carries
# both, and Homebrew installs for your native architecture only. So a plain
# `brew install zstd` on an Apple Silicon Mac gives you an arm64-only libzstd,
# and the universal link will fail with "building for macOS-x86_64 but
# attempting to link with file built for macOS-arm64".
#
# To get a universal build you need a universal libzstd. Either:
#   * build zstd from source twice and `lipo` the results together, or
#   * set ZSTD_LIB to a static universal libzstd.a you already have, or
#   * just run `make` separately on an Intel and an Apple Silicon machine and
#     `lipo -create` the two .dylib files.
#
# Most people should ignore all of this and run `make`.
# ---------------------------------------------------------------------------

UNAME_S  := $(shell uname -s)
UNAME_M  := $(shell uname -m)

CC       ?= cc
PREFIX   ?= /usr/local
DATASET  ?=

BASE_CFLAGS := -O2 -fPIC -Wall -Wextra -std=c99
CFLAGS      ?= $(BASE_CFLAGS)

# --- locate zstd -----------------------------------------------------------
# Precedence: ZSTD_SRC (build from source, link statically) > ZSTD_LIB (an
# explicit library path) > pkg-config > Homebrew > a bare -lzstd.
ifdef ZSTD_SRC
  ZSTD_LOCAL  := $(abspath $(ZSTD_SRC))/lib
  ZSTD_LIB    := $(ZSTD_LOCAL)/libzstd.a
  ZSTD_CFLAGS := -I$(ZSTD_LOCAL)
  ZSTD_LDLIBS := $(ZSTD_LIB)
  ZSTD_STATIC := yes
else
ifndef ZSTD_LIB
  ZSTD_PC := $(shell pkg-config --exists libzstd 2>/dev/null && echo yes)
  ifeq ($(ZSTD_PC),yes)
    ZSTD_CFLAGS := $(shell pkg-config --cflags libzstd)
    ZSTD_LDLIBS := $(shell pkg-config --libs libzstd)
  else
    BREW := $(shell command -v brew 2>/dev/null)
    ifneq ($(BREW),)
      ZSTD_PREFIX := $(shell brew --prefix zstd 2>/dev/null)
    endif
    ifneq ($(ZSTD_PREFIX),)
      ZSTD_CFLAGS := -I$(ZSTD_PREFIX)/include
      ZSTD_LDLIBS := -L$(ZSTD_PREFIX)/lib -lzstd
    else
      ZSTD_CFLAGS :=
      ZSTD_LDLIBS := -lzstd
    endif
  endif
else
  ZSTD_CFLAGS :=
  ZSTD_LDLIBS := $(ZSTD_LIB)
endif
endif

# --- platform --------------------------------------------------------------
ifeq ($(UNAME_S),Darwin)
  LIBNAME    := libtimsdata_tsf.dylib
  COMPATNAME := libtimsdata.dylib
  SHARED     := -dynamiclib -install_name @rpath/$(LIBNAME)
  # macOS ships SQLite in the SDK; no -L needed.
  SQLITE_LDLIBS := -lsqlite3
else
  LIBNAME    := libtimsdata_tsf.so
  COMPATNAME := libtimsdata.so
  SHARED     := -shared -Wl,-soname,$(LIBNAME)
  SQLITE_LDLIBS := -lsqlite3
endif

LDLIBS := $(SQLITE_LDLIBS) $(ZSTD_LDLIBS) -lm
ALL_CFLAGS := $(CFLAGS) $(ZSTD_CFLAGS)

.PHONY: all universal check install compat-name info clean help preflight zstd-static

all: $(LIBNAME)

# Build libzstd.a from a zstd source tree (only when ZSTD_SRC is given).
# zstd's own lib/Makefile does the work; we just pass the compiler and, for a
# universal macOS build, the architecture flags.
ifdef ZSTD_SRC
zstd-static: $(ZSTD_LIB)
$(ZSTD_LIB):
	@test -f "$(ZSTD_LOCAL)/zstd.h" || { \
	  echo "ERROR: $(ZSTD_SRC) does not look like a zstd source tree (no lib/zstd.h)."; \
	  echo "       Download a release from https://github.com/facebook/zstd/releases"; \
	  echo "       and unpack it, then: make ZSTD_SRC=/path/to/zstd-1.5.x"; exit 1; }
	@echo "==> building static libzstd from $(ZSTD_SRC)"
	$(MAKE) -C "$(ZSTD_LOCAL)" libzstd.a CC="$(CC)" CFLAGS="-O3 -fPIC $(ARCHFLAGS)"
preflight: $(ZSTD_LIB)
else
zstd-static:
	@echo "zstd-static needs ZSTD_SRC=/path/to/zstd-source; see the top of this Makefile." >&2; exit 1
endif

# Check the dependencies before compiling, so a missing library produces an
# actionable message rather than a raw linker error.
preflight:
	@printf '#include <zstd.h>\nint main(void){return 0;}\n' > .zstdprobe.c; \
	if ! $(CC) $(ALL_CFLAGS) .zstdprobe.c -o .zstdprobe $(ZSTD_LDLIBS) 2>/dev/null; then \
	  rm -f .zstdprobe .zstdprobe.c; \
	  echo ""; \
	  echo "ERROR: cannot compile and link against zstd."; \
	  echo ""; \
	  if [ "$(UNAME_S)" = "Darwin" ]; then \
	    echo "  Install it with:   brew install zstd"; \
	    echo "  If you have it elsewhere:   make ZSTD_LIB=/path/to/libzstd.a"; \
	  else \
	    echo "  Install it with:   apt install libzstd-dev   (or your package manager)"; \
	  fi; \
	  echo ""; \
	  echo "  'make info' shows what the build is currently detecting."; \
	  echo ""; \
	  exit 1; \
	fi; rm -f .zstdprobe .zstdprobe.c
	@printf '#include <sqlite3.h>\nint main(void){return 0;}\n' > .sqprobe.c; \
	if ! $(CC) $(ALL_CFLAGS) .sqprobe.c -o .sqprobe $(SQLITE_LDLIBS) 2>/dev/null; then \
	  rm -f .sqprobe .sqprobe.c; \
	  echo ""; \
	  echo "ERROR: cannot compile and link against SQLite."; \
	  echo ""; \
	  if [ "$(UNAME_S)" = "Darwin" ]; then \
	    echo "  macOS normally ships SQLite. If this fails, install the Xcode"; \
	    echo "  command line tools:   xcode-select --install"; \
	  else \
	    echo "  Install it with:   apt install libsqlite3-dev"; \
	  fi; \
	  echo ""; \
	  exit 1; \
	fi; rm -f .sqprobe .sqprobe.c

$(LIBNAME): tsf_compat.c tsf_compat.h | preflight
	@echo "==> building $@ for $(UNAME_S)/$(UNAME_M)"
	$(CC) $(ALL_CFLAGS) $(ARCHFLAGS) $(SHARED) -o $@ tsf_compat.c $(LDLIBS)
	@echo
	@echo "Built $@"
	@echo "Load it from your software exactly as you would the vendor library."
	@$(if $(ZSTD_STATIC),echo "zstd is linked statically: this library needs no third-party runtime.",echo "zstd is linked dynamically: the machine running this needs libzstd installed.")
	@echo "Note: m/z is computed from the published calibration equation, so it"
	@echo "agrees with the vendor to about 1 ppm rather than exactly. See README."

universal:
ifneq ($(UNAME_S),Darwin)
	@echo "universal builds are macOS-only; plain 'make' is what you want." >&2
	@exit 1
else
	@echo "==> attempting a universal (arm64 + x86_64) build"
	@echo "    this requires a universal libzstd; see the notes in this Makefile"
	@$(MAKE) --no-print-directory ARCHFLAGS="-arch arm64 -arch x86_64" $(LIBNAME)
	@echo
	@lipo -info $(LIBNAME) 2>/dev/null || true
endif

smoke_test: smoke_test.c $(LIBNAME)
	$(CC) $(ALL_CFLAGS) $(ARCHFLAGS) -o $@ smoke_test.c \
	    -L. -ltimsdata_tsf -Wl,-rpath,$(CURDIR) $(LDLIBS)

check: smoke_test
ifeq ($(DATASET),)
	@echo
	@echo "Built smoke_test. Point it at a dataset to actually check it:"
	@echo "    make check DATASET=/path/to/analysis.d"
	@echo
	@echo "Or, with no data to hand, generate a synthetic one first:"
	@echo "    python ../examples/make_demo_dataset.py /tmp/demo.d"
	@echo "    make check DATASET=/tmp/demo.d"
else
	@echo "==> running smoke test against $(DATASET)"
	@./smoke_test "$(DATASET)"
endif

compat-name: $(LIBNAME)
	@ln -sf $(LIBNAME) $(COMPATNAME)
	@echo "Created symlink $(COMPATNAME) -> $(LIBNAME)."
	@echo
	@echo "This exists so software that dlopen()s the vendor library by its"
	@echo "usual name finds this one instead. It is NOT the vendor library and"
	@echo "is not a drop-in replacement for its numerical output: m/z comes"
	@echo "from the published calibration equation and agrees to about 1 ppm."
	@echo "Use it for local interoperability, not for redistribution under that"
	@echo "name."

install: $(LIBNAME)
	install -d "$(PREFIX)/lib" "$(PREFIX)/include"
	install -m 755 $(LIBNAME) "$(PREFIX)/lib/"
	install -m 644 tsf_compat.h "$(PREFIX)/include/"
	@echo "Installed $(LIBNAME) and tsf_compat.h under $(PREFIX)."
ifeq ($(UNAME_S),Darwin)
	@echo "If your loader cannot find it, set DYLD_LIBRARY_PATH=$(PREFIX)/lib"
	@echo "or link with -Wl,-rpath,$(PREFIX)/lib."
else
	@echo "You may need 'ldconfig' or LD_LIBRARY_PATH=$(PREFIX)/lib."
endif

info:
	@echo "platform     : $(UNAME_S) / $(UNAME_M)"
	@echo "compiler     : $(CC)"
	@echo "library name : $(LIBNAME)"
	@echo "zstd cflags  : $(ZSTD_CFLAGS)"
	@echo "zstd libs    : $(ZSTD_LDLIBS)"
	@echo "zstd linkage : $(if $(ZSTD_STATIC),static (from ZSTD_SRC) -- no runtime dependency,dynamic)"
	@echo "link line    : $(LDLIBS)"
	@echo "arch flags   : $(ARCHFLAGS)"
	@printf "zstd header  : "; \
	  printf '#include <zstd.h>\nint main(void){return 0;}\n' > .zstdprobe.c; \
	  if $(CC) $(ALL_CFLAGS) .zstdprobe.c -o .zstdprobe $(ZSTD_LDLIBS) 2>/dev/null; \
	  then echo "found, links OK"; else echo "NOT USABLE -- install zstd (see top of Makefile)"; fi; \
	  rm -f .zstdprobe .zstdprobe.c
	@printf "sqlite3      : "; \
	  printf '#include <sqlite3.h>\nint main(void){return sqlite3_libversion_number();}\n' > .sqprobe.c; \
	  if $(CC) $(ALL_CFLAGS) .sqprobe.c -o .sqprobe $(SQLITE_LDLIBS) 2>/dev/null; \
	  then echo "found, links OK"; else echo "NOT USABLE -- install SQLite dev headers"; fi; \
	  rm -f .sqprobe .sqprobe.c

help:
	@sed -n '1,40p' Makefile

clean:
	rm -f $(LIBNAME) $(COMPATNAME) smoke_test *.o .zstdprobe* .sqprobe*
