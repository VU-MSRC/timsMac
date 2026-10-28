"""Reader for the ``analysis.tsf_bin`` binary spectrum file.

The structure of this file is documented by `timsrust-tsf`_ (MannLabs, MIT),
on which this reader is based; see LICENSE for its terms. Behaviour not covered
there was established by testing against the vendor library through its
documented public C API, and is expressed here in the implementation rather
than restated as a specification. See PROVENANCE.md.

Spectra are located by ``Frames.TimsId`` and decompressed with zstd. Line
(centroided) spectra yield index, intensity and -- where the acquisition
recorded them -- peak width arrays; profile spectra yield a dense array in
which position is the index.

.. _timsrust-tsf: https://github.com/MannLabs/timsrust
"""

from __future__ import annotations

import contextlib
import struct

import numpy as np
import zstandard as zstd

__all__ = [
    "MAX_SUBBLOCK_SIZE",
    "BinaryFile",
    "TsfBinError",
    "decode_line_block",
    "decode_profile_block",
    "encode_line_block",
    "encode_profile_block",
    "encode_spectrum_block",
    "transpose_to_byte_planes",
]

#: Largest sub-block this reader will attempt to read (1 GiB). A declared
#: size beyond this is treated as corruption rather than allocated.
MAX_SUBBLOCK_SIZE = 0x40000000

_HEADER = struct.Struct("<II")
_DCTX = zstd.ZstdDecompressor()


class TsfBinError(Exception):
    """Raised when a .tsf_bin block cannot be parsed."""


def _decompress(buf: bytes, what: str) -> bytes:
    """Decompress one raw zstd frame.

    Frames written by the acquisition software declare their content size; we
    try that path first and fall back to streaming for any that do not.
    """
    try:
        return _DCTX.decompress(buf)
    except zstd.ZstdError:
        pass
    try:
        dobj = zstd.ZstdDecompressor().decompressobj()
        return dobj.decompress(buf) + dobj.flush()
    except zstd.ZstdError as exc:
        raise TsfBinError(f"error decompressing {what} data: {exc}") from exc


# -- decoding -------------------------------------------------------------


def decode_line_block(raw: bytes, num_peaks: int, with_width: bool | None = None):
    """Decode an already-decompressed line sub-block.

    Returns ``(index, intensity, width_or_None)``.
    """
    n = int(num_peaks)
    if n < 0:
        raise TsfBinError("NumPeaks is NULL/negative; no line spectrum stored")
    if n == 0:
        return _empty_line(bool(with_width))

    if with_width is None:
        with_width = len(raw) >= n * 16
    need = n * (16 if with_width else 12)
    if len(raw) < need:
        raise TsfBinError(f"line sub-block too short: {len(raw)} bytes, need {need} for {n} peaks")

    index = np.frombuffer(raw, dtype="<f8", count=n, offset=0).copy()
    inten = np.frombuffer(raw, dtype="<f4", count=n, offset=n * 8).copy()
    width = np.frombuffer(raw, dtype="<f4", count=n, offset=n * 12).copy() if with_width else None
    return index, inten, width


def _empty_line(with_width: bool):
    return (
        np.empty(0, dtype="<f8"),
        np.empty(0, dtype="<f4"),
        np.empty(0, dtype="<f4") if with_width else None,
    )


def decode_profile_block(raw: bytes) -> np.ndarray:
    """Decode an already-decompressed profile sub-block (byte planes)."""
    n = len(raw) // 4  # a non-multiple-of-4 tail is ignored
    if n == 0:
        return np.empty(0, dtype=np.uint32)
    planes = np.frombuffer(raw, dtype=np.uint8, count=4 * n).reshape(4, n)
    out = planes[0].astype(np.uint32)
    out |= planes[1].astype(np.uint32) << np.uint32(8)
    out |= planes[2].astype(np.uint32) << np.uint32(16)
    out |= planes[3].astype(np.uint32) << np.uint32(24)
    return out


# -- file access ----------------------------------------------------------


class BinaryFile:
    """An open ``analysis.tsf_bin``.

    Holds one file handle for the lifetime of the object, which matters for
    MALDI imaging analyses where a single file can contain >10^5 spectra.
    """

    def __init__(self, path: str):
        self.path = str(path)
        self._f = open(self.path, "rb")

    # -- lifecycle --
    def close(self) -> None:
        if getattr(self, "_f", None) is not None:
            self._f.close()
            self._f = None  # type: ignore[assignment]

    def __enter__(self) -> BinaryFile:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.close()

    # -- header --
    def block_header(self, tims_id: int | None) -> tuple[int, int, int]:
        """Return ``(total_size, line_size, profile_size)`` for one block."""
        if tims_id is None or tims_id < 0:
            return (0, 0, 0)
        self._f.seek(tims_id)
        hdr = self._f.read(8)
        if len(hdr) != 8:
            raise TsfBinError(f"error reading binary header at offset {tims_id}")
        total, line_size = _HEADER.unpack(hdr)
        # [rust blobs.rs:91] checks chunk_padded >= 8 and >= compressed_len.
        # The stricter relation total >= 8 + line_size must hold, because
        # total delimits the profile sub-block.
        if total < 8 or line_size > total - 8:
            raise TsfBinError(f"corrupt block header at offset {tims_id}: " f"total={total} line={line_size}")
        return (total, line_size, total - 8 - line_size)

    def _sub_block(self, tims_id: int, which: str) -> bytes:
        total, line_size, prof_size = self.block_header(tims_id)
        if total == 0:
            return b""
        if which == "line":
            offset, size = tims_id + 8, line_size
        else:
            offset, size = tims_id + 8 + line_size, prof_size
        if size == 0:
            return b""
        if size > MAX_SUBBLOCK_SIZE:
            raise TsfBinError(f"compressed {which} spectrum size too large ({size})")
        self._f.seek(offset)
        buf = self._f.read(size)
        if len(buf) != size:
            raise TsfBinError(
                f"error reading {which} frame data at offset {offset}: " f"wanted {size} bytes, got {len(buf)}"
            )
        return buf

    # -- reads --
    def read_line_spectrum(
        self,
        tims_id: int | None,
        num_peaks: int | None,
        want_width: bool | None = None,
    ):
        """Read one line spectrum as ``(index, intensity, width_or_None)``.

        ``want_width`` defaults to inferring from the payload length rather
        than trusting ``GlobalMetadata.HasLineSpectraPeakWidth``, so that a
        flag inconsistent with the payload cannot cause an out-of-bounds read.
        """
        if tims_id is None or tims_id < 0 or num_peaks is None or int(num_peaks) <= 0:
            return _empty_line(bool(want_width))
        raw = self._sub_block(tims_id, "line")
        if not raw:
            return _empty_line(bool(want_width))
        return decode_line_block(_decompress(raw, "line"), num_peaks, want_width)

    def read_profile_spectrum(self, tims_id: int | None) -> np.ndarray:
        """Read one profile spectrum as a dense uint32 array (position == index)."""
        if tims_id is None or tims_id < 0:
            return np.empty(0, dtype=np.uint32)
        raw = self._sub_block(tims_id, "profile")
        if not raw:
            return np.empty(0, dtype=np.uint32)
        return decode_profile_block(_decompress(raw, "profile"))


# -- encoding (used by the test suite to build synthetic fixtures) --------


def transpose_to_byte_planes(values) -> bytes:
    """uint32 sequence -> the byte-plane layout of a profile sub-block.

    Inverse of :func:`decode_profile_block`.
    """
    v = np.asarray(values, dtype="<u4")
    if v.size == 0:
        return b""
    return np.frombuffer(v.tobytes(), dtype=np.uint8).reshape(-1, 4).T.tobytes()


def encode_line_block(index, intensity, width=None, level: int = 3) -> bytes:
    if len(index) != len(intensity):
        raise ValueError("index/intensity length mismatch")
    if width is not None and len(width) != len(index):
        raise ValueError("width length mismatch")
    parts = [
        np.asarray(index, dtype="<f8").tobytes(),
        np.asarray(intensity, dtype="<f4").tobytes(),
    ]
    if width is not None:
        parts.append(np.asarray(width, dtype="<f4").tobytes())
    return zstd.ZstdCompressor(level=level).compress(b"".join(parts))


def encode_profile_block(values, level: int = 3) -> bytes:
    return zstd.ZstdCompressor(level=level).compress(transpose_to_byte_planes(values))


def encode_spectrum_block(index=(), intensity=(), width=None, profile=None, level: int = 3) -> bytes:
    """Build one complete per-spectrum block (8-byte header + both sub-blocks).

    ``profile=None`` omits the profile sub-block entirely; ``profile=[]``
    stores an empty zstd frame, which is how "no profile" is represented.
    """
    line = encode_line_block(index, intensity, width, level) if len(index) else b""
    prof = b"" if profile is None else encode_profile_block(profile, level)
    return _HEADER.pack(8 + len(line) + len(prof), len(line)) + line + prof
