"""Ion-mobility calibration for TDF data: scan number -> 1/K0 -> CCS.

Provenance
----------
Nothing here was obtained from the vendor library, which records its own
mobility calibration in a `TimsCalibration` table whose coefficient layout is
not published and which this package does not read. Two calibrations are
offered, both built from published material and the file's own contents:

* :meth:`MobilityCalibration.from_reference_ions` -- **external calibration**:
  the standard practice in ion-mobility mass spectrometry. Ions of known 1/K0
  (for example the Agilent ESI tune mix, whose values are tabulated in the CCS
  compendium of Stow *et al.*, *Anal. Chem.* 2017, 89, 9048) are located in
  the data and a low-order polynomial in scan number is fitted through them.
  The physical basis is that TIMS elutes ions as the trapping field is ramped
  linearly in time, and 1/K0 is to first order linear in that field (Michelmann
  *et al.*, *J. Am. Soc. Mass Spectrom.* 2015, 26, 14; Ridgeway *et al.*,
  *Int. J. Mass Spectrom.* 2018, 425, 22). Measured against the vendor library
  on a direct-infusion tune-mix run: median 0.25 % inside the calibrant span,
  against an instrument-stated calibration SD of 0.26 %.

  When the file's own calibrant (recorded in the documented `CalibrationInfo`
  table as `ReferencePeakMasses` / `ReferencePeakMobilities`) is present in the
  spectra, :func:`find_reference_ions` locates it automatically.

* :meth:`MobilityCalibration.from_acquisition_range` -- **approximation**:
  linear interpolation between the documented `OneOverK0AcqRangeLower/Upper`
  metadata across the scan range. This is what `timsrust` ships. It needs no
  ions in the data, but the real ramp overshoots the nominal range, so expect
  errors of 2-4 % (median 2.1 % and 3.5 % on two files). It is provided so
  that files without usable calibrant ions can still be read, and it warns.

CCS is obtained from 1/K0 with the Mason-Schamp equation, textbook physics.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    "MobilityCalibration",
    "ReferenceIonFit",
    "ccs_to_oneoverk0",
    "find_reference_ions",
    "oneoverk0_to_ccs",
]

# Mason-Schamp: CCS = (3 z e / 16 N0) * sqrt(2 pi / (mu kB T)) / K0, with K0 the
# reduced mobility at standard number density N0 (Loschmidt constant). The
# prefactor is computed here from CODATA 2018 values so that it is traceable;
# it evaluates to 18509.87, the figure usually quoted in the TIMS literature.
# Buffer gas is N2. The TIMS community evaluates the formula at T = 305 K, and
# published CCS tables (e.g. Stow et al. 2017) are on that convention.
_E = 1.602176634e-19  # C
_KB = 1.380649e-23  # J/K
_AMU = 1.66053906660e-27  # kg
_N0 = 2.686780111e25  # 1/m^3 at 273.15 K, 101.325 kPa
_N2_MASS = 28.0134  # Da
_MS_CONST = 1e24 * 3.0 * _E / (16.0 * _N0) * (2.0 * np.pi / (_AMU * _KB)) ** 0.5


def _reduced_mass(mz: np.ndarray, charge: int) -> np.ndarray:
    m_ion = mz * charge
    return m_ion * _N2_MASS / (m_ion + _N2_MASS)


def oneoverk0_to_ccs(oneoverk0, mz, charge: int = 1, temperature_k: float = 305.0):
    """Mason-Schamp: CCS (Å²) from 1/K0 (V·s/cm²), m/z and charge."""
    oneoverk0 = np.asarray(oneoverk0, dtype=np.float64)
    mz = np.asarray(mz, dtype=np.float64)
    mu = _reduced_mass(mz, charge)
    return _MS_CONST * charge * oneoverk0 / np.sqrt(mu * temperature_k)


def ccs_to_oneoverk0(ccs, mz, charge: int = 1, temperature_k: float = 305.0):
    ccs = np.asarray(ccs, dtype=np.float64)
    mz = np.asarray(mz, dtype=np.float64)
    mu = _reduced_mass(mz, charge)
    return ccs * np.sqrt(mu * temperature_k) / (_MS_CONST * charge)


@dataclass
class ReferenceIonFit:
    """What :func:`find_reference_ions` found for one reference ion."""

    name: str
    mz: float
    oneoverk0: float
    apex_scan: float
    fwhm_scans: float
    observed_mz: float
    intensity: float
    found: bool
    snr: float = float("nan")  # apex height over the profile's median level

    @property
    def mz_error_ppm(self) -> float:
        return (self.observed_mz - self.mz) / self.mz * 1e6


@dataclass
class MobilityCalibration:
    """scan number <-> 1/K0 for one analysis (one ramp)."""

    coefficients: np.ndarray  # polynomial in scan, highest power first
    num_scans: int
    method: str  # "reference_ions" | "acquisition_range"
    scan_min: float  # span over which the fit is anchored
    scan_max: float
    fit_residual_percent: float = 0.0
    reference_ions: list = field(default_factory=list)
    warn_outside_range: bool = True
    _warned: bool = field(default=False, repr=False, compare=False)

    # -- constructors -----------------------------------------------------

    @classmethod
    def from_reference_ions(
        cls,
        apex_scans: Sequence[float],
        oneoverk0: Sequence[float],
        num_scans: int,
        order: int = 1,
        warn_outside_range: bool = True,
        reference_ions: list | None = None,
    ) -> MobilityCalibration:
        """External calibration from ions of known 1/K0 observed at known scans.

        ``order=1`` (default) is the physically motivated linear relation;
        ``order=2`` fits a quadratic and generalises less safely beyond the
        calibrant span. At least ``order + 1`` ions are required.
        """
        s = np.asarray(apex_scans, dtype=np.float64)
        k = np.asarray(oneoverk0, dtype=np.float64)
        if s.size < order + 1:
            raise ValueError(f"need at least {order + 1} reference ions for order {order}")
        coef = np.polyfit(s, k, order)
        resid = np.polyval(coef, s) - k
        return cls(
            coef,
            int(num_scans),
            "reference_ions",
            float(s.min()),
            float(s.max()),
            float(np.max(np.abs(resid)) / k.mean() * 100.0),
            list(reference_ions or []),
            warn_outside_range,
        )

    @classmethod
    def from_acquisition_range(
        cls,
        oneoverk0_lower: float,
        oneoverk0_upper: float,
        num_scans: int,
        warn_outside_range: bool = True,
        warn: bool = True,
    ) -> MobilityCalibration:
        """Approximation: linear between the documented acquisition-range
        endpoints, high 1/K0 at scan 0. Expect 2-4 % error; see module docs.
        ``warn=False`` silences the notice (TdfData issues its own).
        """
        n = int(num_scans)
        slope = (oneoverk0_lower - oneoverk0_upper) / max(n - 1, 1)
        coef = np.array([slope, oneoverk0_upper])
        if warn:
            warnings.warn(
                "Mobility calibrated from the nominal acquisition range only "
                "(no reference ions): expect 1/K0 errors of a few percent. Pass "
                "reference ions, or acquire with calibrant present, for ~0.3 %.",
                RuntimeWarning,
                stacklevel=3,
            )
        return cls(
            coef,
            n,
            "acquisition_range",
            0.0,
            float(n - 1),
            float("nan"),
            [],
            warn_outside_range,
        )

    # -- conversion -------------------------------------------------------

    def scan_to_oneoverk0(self, scans):
        scalar = np.isscalar(scans)
        s = np.atleast_1d(np.asarray(scans, dtype=np.float64))
        k = np.polyval(self.coefficients, s)
        if self.warn_outside_range and not self._warned and self.method == "reference_ions" and s.size:
            out = (s < self.scan_min - 0.5) | (s > self.scan_max + 0.5)
            if out.any():
                self._warned = True  # once per calibration, not once per call
                warnings.warn(
                    f"scans outside the span of the reference ions "
                    f"({self.scan_min:.0f}-{self.scan_max:.0f}) were converted; the "
                    "mobility calibration extrapolates there and accuracy degrades. "
                    "(Reported once per calibration.)",
                    RuntimeWarning,
                    stacklevel=2,
                )
        return float(k[0]) if scalar else k

    def oneoverk0_to_scan(self, oneoverk0):
        scalar = np.isscalar(oneoverk0)
        k = np.atleast_1d(np.asarray(oneoverk0, dtype=np.float64))
        c = self.coefficients
        if c.size == 2:
            s = (k - c[1]) / c[0]
        else:
            # solve the quadratic/cubic numerically, seeded from the linear part
            s = np.empty_like(k)
            for i, kk in enumerate(k):
                roots = np.roots(np.r_[c[:-1], c[-1] - kk])
                real = roots[np.isreal(roots)].real
                inside = real[(real >= -1) & (real <= self.num_scans)]
                s[i] = inside[np.argmin(np.abs(inside - self.num_scans / 2))] if inside.size else np.nan
        return float(s[0]) if scalar else s

    def __repr__(self) -> str:
        if self.method == "reference_ions":
            return (
                f"MobilityCalibration(reference_ions={len(self.reference_ions) or '?'}, "
                f"order={self.coefficients.size - 1}, scans {self.scan_min:.0f}-{self.scan_max:.0f}, "
                f"fit_residual={self.fit_residual_percent:.3f}%)"
            )
        return (
            f"MobilityCalibration(acquisition_range approximation, "
            f"{np.polyval(self.coefficients, 0):.3f} -> "
            f"{np.polyval(self.coefficients, self.num_scans - 1):.3f})"
        )


# -- locating reference ions in the data -----------------------------------


def find_reference_ions(
    tdf,
    mz_list: Sequence[float],
    oneoverk0_list: Sequence[float],
    names: Sequence[str] | None = None,
    ppm: float = 15.0,
    frame_ids: Sequence[int] | None = None,
    min_intensity: float = 0.0,
) -> list[ReferenceIonFit]:
    """Locate ions of known m/z and 1/K0 in the data and return their scan apex.

    For each ion, all peaks within ``ppm`` of its m/z (on this package's m/z
    axis) are accumulated into an intensity-vs-scan profile over ``frame_ids``
    (default: all MS1 frames); the apex is the intensity-weighted centroid of
    the region above half maximum. Ions with no signal are returned with
    ``found=False``.
    """
    names = list(names) if names is not None else [f"{m:.4f}" for m in mz_list]
    mz_list = [float(m) for m in mz_list]
    frames = list(frame_ids) if frame_ids is not None else [f.id for f in tdf.frames if f.msms_type == 0]
    ns = tdf.meta.num_scans
    cal = tdf.calibration

    lo = cal.mz_to_index(np.array(mz_list) * (1 - ppm * 1e-6))
    hi = cal.mz_to_index(np.array(mz_list) * (1 + ppm * 1e-6))
    acc = np.zeros((len(mz_list), ns))
    tof_w = np.zeros(len(mz_list))
    tof_i = np.zeros(len(mz_list))
    for fid in frames:
        fr = tdf.read_frame(fid)
        if fr.num_peaks == 0:
            continue
        tof = fr.tof_indices.astype(np.float64)
        sc = fr.scan
        it = fr.intensities.astype(np.float64)
        for j in range(len(mz_list)):
            m = (tof >= lo[j]) & (tof <= hi[j])
            if m.any():
                np.add.at(acc[j], sc[m], it[m])
                tof_w[j] += (tof[m] * it[m]).sum()
                tof_i[j] += it[m].sum()

    out = []
    for j, (nm, mz, k0) in enumerate(zip(names, mz_list, oneoverk0_list)):
        a = acc[j]
        if a.sum() <= min_intensity or a.max() <= 0:
            out.append(ReferenceIonFit(nm, mz, float(k0), np.nan, np.nan, np.nan, float(a.sum()), False))
            continue
        pk = int(np.argmax(a))
        half = a[pk] / 2.0
        l = pk
        r = pk
        while l > 0 and a[l - 1] >= half:
            l -= 1
        while r < ns - 1 and a[r + 1] >= half:
            r += 1
        seg = slice(max(l - 3, 0), min(r + 4, ns))
        s = np.arange(ns, dtype=np.float64)[seg]
        w = a[seg]
        apex = float((s * w).sum() / w.sum())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obs_mz = float(cal.index_to_mz(tof_w[j] / tof_i[j])) if tof_i[j] else np.nan
        # apex height over the background: the median of the non-empty scans
        # outside the peak (peak = +/- 3 FWHM, at least +/- 10 scans). A real
        # calibrant stands far above it; a noise pick does not. With no
        # background at all the ratio is infinite.
        margin = max(10, 3 * (r - l + 1))
        bg = np.concatenate([a[: max(l - margin, 0)], a[min(r + margin + 1, ns) :]])
        bg = bg[bg > 0]
        snr = float(a[pk] / np.median(bg)) if bg.size else float("inf")
        out.append(
            ReferenceIonFit(
                nm,
                mz,
                float(k0),
                apex,
                float(r - l + 1),
                obs_mz,
                float(a.sum()),
                True,
                snr,
            )
        )
    return out
