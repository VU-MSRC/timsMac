"""Reader for the ``analysis.tdf_bin`` binary frame file (TDF: ion-mobility data).

Provenance
----------
The layout implemented here is the one published by `timsrust`_ (MannLabs) in
``crates/timsrust-tdf/src/frame_reader/compression2.rs`` and, independently, by
OpenTIMS. Nothing here was obtained from the vendor library; see PROVENANCE.md.

Layout (``GlobalMetadata.TimsCompressionType == 2``)
-----------------------------------------------------
At byte offset ``Frames.TimsId``::

    +0   uint32 LE   byte_count   whole block, including this 8-byte header
    +4   uint32 LE   (scan count; the same value is repeated inside the payload)
    +8   byte_count - 8 bytes of zstd-compressed payload

The decompressed payload is a byte-plane-transposed uint32 array: with
``n = len // 4``, value ``i`` is ``p0[i] | p1[i]<<8 | p2[i]<<16 | p3[i]<<24``
where ``pk = bytes[k*n:(k+1)*n]``. In that array:

    u[0]                    scan_count
    u[1 .. scan_count)      size of scans 0 .. scan_count-2, in uint32 units
                            (peaks in that scan = size / 2); the last scan
                            takes whatever remains
    u[scan_count ..]        interleaved (tof_delta, intensity) pairs, scan by
                            scan; within a scan the TOF index is the running
                            sum of the deltas, minus one

``peak_count = (n - scan_count) / 2``.

.. _timsrust: https://github.com/MannLabs/timsrust
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

import numpy as np
import zstandard as zstd

__all__ = [
    "Frame",
    "TdfBinError",
    "TdfBinaryFile",
    "decode_frame_payload",
    "encode_frame_payload",
]

_HEADER = struct.Struct("<II")
_DCTX = zstd.ZstdDecompressor()

#: Largest block this reader will attempt (1 GiB).
MAX_BLOCK_SIZE = 0x40000000


class TdfBinError(Exception):
    """Raised when a .tdf_bin block cannot be parsed."""


@dataclass
class Frame:
    """One decoded frame: parallel arrays of (scan, TOF index, intensity).

    ``scan_offsets`` has ``num_scans + 1`` entries; the peaks of scan ``s`` are
    ``tof_indices[scan_offsets[s]:scan_offsets[s+1]]``. ``scan`` is the same
    information expanded to one entry per peak, for convenience.
    """

    num_scans: int
    scan_offsets: np.ndarray  # int64, len num_scans + 1
    tof_indices: np.ndarray  # uint32
    intensities: np.ndarray  # uint32

    @property
    def scan(self) -> np.ndarray:
        counts = np.diff(self.scan_offsets)
        return np.repeat(np.arange(self.num_scans, dtype=np.int32), counts)

    @property
    def num_peaks(self) -> int:
        return int(self.tof_indices.size)


def _unshuffle(raw: bytes) -> np.ndarray:
    n = len(raw) // 4
    if n == 0 or len(raw) % 4:
        raise TdfBinError(f"payload length {len(raw)} is not a multiple of 4")
    planes = np.frombuffer(raw, dtype=np.uint8, count=4 * n).reshape(4, n)
    out = planes[0].astype(np.uint32)
    out |= planes[1].astype(np.uint32) << np.uint32(8)
    out |= planes[2].astype(np.uint32) << np.uint32(16)
    out |= planes[3].astype(np.uint32) << np.uint32(24)
    return out


def decode_frame_payload(raw: bytes) -> Frame:
    """Decode an already-decompressed frame payload."""
    u = _unshuffle(raw)
    scan_count = int(u[0])
    if scan_count < 1 or scan_count > u.size:
        raise TdfBinError(f"implausible scan count {scan_count}")
    peak_count = (u.size - scan_count) // 2
    if (u.size - scan_count) % 2:
        raise TdfBinError("payload does not hold whole (tof, intensity) pairs")

    # per-scan peak counts: scans 0..scan_count-2 are given, the last is the remainder
    sizes = (u[1:scan_count] // 2).astype(np.int64)
    scan_offsets = np.empty(scan_count + 1, dtype=np.int64)
    scan_offsets[0] = 0
    np.cumsum(sizes, out=scan_offsets[1:scan_count])
    scan_offsets[scan_count] = peak_count
    if scan_offsets[scan_count - 1] > peak_count:
        raise TdfBinError("scan sizes exceed the peak count")

    pairs = u[scan_count : scan_count + 2 * peak_count].reshape(peak_count, 2)
    tof_deltas = pairs[:, 0].astype(np.int64)
    intensities = pairs[:, 1].copy()

    # Running sum of deltas within each scan, minus one. Take the global
    # running sum, then subtract from every peak the total accumulated before
    # its scan began -- correct for empty scans too, since those contribute
    # nothing to the global sum.
    tof = np.cumsum(tof_deltas)
    counts = np.diff(scan_offsets)
    starts = scan_offsets[:-1]
    before = np.where(starts > 0, tof[np.maximum(starts - 1, 0)], 0)
    tof = tof - np.repeat(before, counts)
    tof_indices = (tof - 1).astype(np.uint32)

    return Frame(scan_count, scan_offsets, tof_indices, intensities)


class TdfBinaryFile:
    """An open ``analysis.tdf_bin``; holds one file handle."""

    def __init__(self, path: str):
        self.path = str(path)
        self._f = open(self.path, "rb")

    def close(self) -> None:
        if getattr(self, "_f", None) is not None:
            self._f.close()
            self._f = None  # type: ignore[assignment]

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def read_frame(self, tims_id: int | None) -> Frame:
        if tims_id is None or tims_id < 0:
            return Frame(0, np.zeros(1, np.int64), np.empty(0, np.uint32), np.empty(0, np.uint32))
        self._f.seek(tims_id)
        hdr = self._f.read(8)
        if len(hdr) != 8:
            raise TdfBinError(f"error reading frame header at offset {tims_id}")
        byte_count, _hdr_scans = _HEADER.unpack(hdr)
        if byte_count < 8:
            raise TdfBinError(f"corrupt frame header at offset {tims_id}: byte_count={byte_count}")
        if byte_count > MAX_BLOCK_SIZE:
            raise TdfBinError(f"frame block too large ({byte_count})")
        comp = self._f.read(byte_count - 8)
        if len(comp) != byte_count - 8:
            raise TdfBinError(f"short read of frame data at offset {tims_id}")
        if not comp:
            return Frame(
                int(_hdr_scans),
                np.zeros(int(_hdr_scans) + 1, np.int64),
                np.empty(0, np.uint32),
                np.empty(0, np.uint32),
            )
        try:
            raw = _DCTX.decompress(comp)
        except zstd.ZstdError as exc:
            raise TdfBinError(f"error decompressing frame at offset {tims_id}: {exc}") from exc
        return decode_frame_payload(raw)


# -- encoder, for tests ---------------------------------------------------


def encode_frame_payload(num_scans: int, scan_of_peak, tof_indices, intensities, level: int = 3) -> bytes:
    """Build one complete frame block (header + zstd payload) from peaks.

    Peaks must be sorted by scan, then by TOF index within a scan.
    """
    scan_of_peak = np.asarray(scan_of_peak, dtype=np.int64)
    tof_indices = np.asarray(tof_indices, dtype=np.int64)
    intensities = np.asarray(intensities, dtype=np.uint32)
    counts = np.bincount(scan_of_peak, minlength=num_scans).astype(np.int64)
    # deltas: first peak in a scan is tof+1, subsequent are differences
    deltas = np.empty_like(tof_indices)
    pos = 0
    for s in range(num_scans):
        c = counts[s]
        if c:
            seg = tof_indices[pos : pos + c]
            deltas[pos] = seg[0] + 1
            deltas[pos + 1 : pos + c] = np.diff(seg)
            pos += c
    u = np.empty(num_scans + 2 * tof_indices.size, dtype=np.uint32)
    u[0] = num_scans
    u[1:num_scans] = (counts[:-1] * 2).astype(np.uint32)
    u[num_scans::2] = deltas.astype(np.uint32)
    u[num_scans + 1 :: 2] = intensities
    shuffled = np.frombuffer(u.tobytes(), dtype=np.uint8).reshape(-1, 4).T.tobytes()
    comp = zstd.ZstdCompressor(level=level).compress(shuffled)
    return _HEADER.pack(8 + len(comp), num_scans) + comp
