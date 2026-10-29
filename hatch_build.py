"""Optionally bundle the standalone TSF C library in platform wheels."""

from __future__ import annotations

import os
import shutil
import subprocess  # ruff: ignore[suspicious-subprocess-import]
import sys
import sysconfig
from pathlib import Path
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from packaging.tags import mac_platforms


class CapiBuildHook(BuildHookInterface):
    """Compile the C API and add its library and header to native wheels."""

    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        """Run the existing Makefile only when a native wheel is requested."""
        if os.environ.get("PYTSF_BUILD_CAPI") != "1":
            return
        if sys.platform not in {"darwin", "linux"}:
            raise RuntimeError("The C API Makefile supports Linux and macOS only.")

        source = Path(self.root) / "capi"
        temporary = Path(self.root) / "build" / "capi"
        temporary.mkdir(parents=True, exist_ok=True)
        for name in ("Makefile", "tsf_compat.c", "tsf_compat.h"):
            shutil.copy2(source / name, temporary / name)
        subprocess.run(["make", "-B", "-C", str(temporary)], check=True)  # ruff: ignore[subprocess-without-shell-equals-true, start-process-with-partial-path]

        suffix = ".dylib" if sys.platform == "darwin" else ".so"
        library = f"libtimsdata_tsf{suffix}"
        build_data["force_include"][str(temporary / library)] = f"pytsf/{library}"
        build_data["force_include"][str(source / "tsf_compat.h")] = "pytsf/tsf_compat.h"
        build_data["pure_python"] = False
        platform_tag = sysconfig.get_platform().replace("-", "_").replace(".", "_")
        if sys.platform == "darwin":
            target = os.environ.get("MACOSX_DEPLOYMENT_TARGET")
            deployment = None
            if target:
                major, _, minor = target.partition(".")
                deployment = (int(major), int(minor or "0"))
            platform_tag = next(mac_platforms(version=deployment))
        build_data["tag"] = f"py3-none-{platform_tag}"
