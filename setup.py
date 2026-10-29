"""Optionally bundle the standalone TSF C library in platform wheels."""

from __future__ import annotations

import os
import shutil
import subprocess  # ruff: ignore[suspicious-subprocess-import]
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_ext import build_ext

BUILD_CAPI = os.environ.get("PYTSF_BUILD_CAPI") == "1"
LIBRARY_SUFFIX = ".dylib" if sys.platform == "darwin" else ".so"


class BuildCapi(build_ext):
    """Build the C ABI library with its existing Makefile."""

    def get_ext_filename(self, fullname: str) -> str:
        """Use the standalone library name without a Python ABI suffix."""
        return str(Path(*fullname.split(".")).with_suffix(LIBRARY_SUFFIX))

    def build_extension(self, ext: Extension) -> None:
        """Build in a temporary directory and install the library and header."""
        if sys.platform not in {"darwin", "linux"}:
            raise RuntimeError("The C API Makefile supports Linux and macOS only.")

        source = Path(__file__).parent / "capi"
        temporary = Path(self.build_temp).resolve() / "capi"
        temporary.mkdir(parents=True, exist_ok=True)
        for name in ("Makefile", "tsf_compat.c", "tsf_compat.h"):
            shutil.copy2(source / name, temporary / name)
        subprocess.run(["make", "-B", "-C", str(temporary)], check=True)  # ruff: ignore[subprocess-without-shell-equals-true, start-process-with-partial-path]

        destination = Path(self.get_ext_fullpath(ext.name))
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(temporary / f"libtimsdata_tsf{LIBRARY_SUFFIX}", destination)
        shutil.copy2(source / "tsf_compat.h", destination.parent / "tsf_compat.h")


class CapiWheel(bdist_wheel):
    """Tag the standalone C library as independent of the Python ABI."""

    def get_tag(self) -> tuple[str, str, str]:
        """Keep the platform tag while allowing any supported Python 3."""
        python, abi, platform = super().get_tag()
        return ("py3", "none", platform) if BUILD_CAPI else (python, abi, platform)


setup(
    ext_modules=[Extension("pytsf.libtimsdata_tsf", sources=["capi/tsf_compat.c"])] if BUILD_CAPI else [],
    cmdclass={"build_ext": BuildCapi, "bdist_wheel": CapiWheel},
)
