"""Independent m/z calibration: the published TOF equation fitted to the
calibration points the instrument recorded in the file.

Provenance
----------
This module is built entirely from published and vendor-documented material:

* **The equation.** The three-point time-of-flight calibration

      t = ml2 + sqrt(1e12/ml1 * m) + ml3 * m

  as implemented by OpenMS's ``TOFCalibration``, whose documentation states
  that the procedure "is very similar to the one described in" Gobom, J.,
  Mueller, M., Egelhofer, V., Theiss, D., Lehrach, H. & Nordhoff, E., "A
  calibration method that simplifies and improves accurate determination of
  peptide molecular masses by MALDI-TOF MS", *Analytical Chemistry* 2002, 74
  (15), 3915-3923. ``order=2`` drops the ``ml3`` term, giving the
  classical two-point form; ``order=4`` adds an ``ml4 * m**1.5`` term.

* **The calibration points.** ``CalibrationInfo`` -- a table Bruker documents
  in the schema reference distributed with the SDK as the "List of calibration
  information for this analysis" -- stores ``ReferencePeakMasses`` (the known
  masses of the calibrant used) and ``MeasuredTimesOfFlight`` (the flight times
  at which those peaks were observed). Those are exactly the (t, m) pairs the
  equation above is fitted to.

* **The time axis.** ``MzCalibration.DigitizerTimebase`` and
  ``DigitizerDelay``, documented columns whose names carry their ordinary
  digitiser meaning: flight time is affine in sample number, so
  ``t = delay + index * timebase``.

Nothing here was obtained by examining the vendor library, and ``pytsf``
deliberately does not reproduce the vendor's own calibration algorithm, which
its documentation describes as proprietary. See PROVENANCE.md.

Accuracy
--------
This is an *external calibration*: one fit per analysis, accurate across the
mass range the calibrant covered and degrading when extrapolating beyond it,
exactly as any external calibration does. ``index_to_mz`` can warn when asked
to extrapolate (``warn_outside_range``). Measured figures against the vendor
library on real data are in the README.
"""

from __future__ import annotations

import sqlite3
import struct
import warnings

import numpy as np

from .meta import TsfError

__all__ = ["FittedCalibration", "CalibrationPointsUnavailable"]

SCALE = 1e12


class CalibrationPointsUnavailable(TsfError):
    """The analysis records no reference masses / measured flight times."""


def _doubles(b) -> np.ndarray:
    if not isinstance(b, (bytes, bytearray, memoryview)):
        raise CalibrationPointsUnavailable(
            "Calibration point data is not a binary array.")
    b = bytes(b)
    if len(b) < 8 or len(b) % 8:
        raise CalibrationPointsUnavailable(
            f"Calibration point blob has unusable length {len(b)}.")
    return np.array(struct.unpack(f"<{len(b) // 8}d", b))


class FittedCalibration:
    """Fit the published TOF equation to a file's own calibration points.

    Parameters
    ----------
    order
        3 (default) for the published three-point equation, 2 for the
        two-point form, 4 to add one further term. More terms fit the
        calibrant better but generalise less safely beyond it.
    warn_outside_range
        Emit a warning when converting indices whose m/z falls outside the
        span of the calibrant masses, where accuracy degrades.
    """

    def __init__(self, sqlite_path: str, order: int = 3,
                 warn_outside_range: bool = True):
        if order not in (2, 3, 4):
            raise ValueError("order must be 2, 3 or 4")
        self.order = order
        self.warn_outside_range = warn_outside_range
        self._warned = False
        self._axis_cache: dict[int, np.ndarray] = {}

        uri = ("file:" + str(sqlite_path).replace("?", "%3f").replace("#", "%23")
               + "?mode=ro")
        con = sqlite3.connect(uri, uri=True)
        try:
            self._load(con)
        finally:
            con.close()
        self._fit()

    # -- loading ----------------------------------------------------------

    def _load(self, con: sqlite3.Connection) -> None:
        # One fit covers the whole analysis, so the analysis must be
        # single-polarity and reference a single calibration row. Refuse
        # anything else rather than silently applying one calibration to all.
        polarities = [r[0] for r in con.execute(
            "SELECT DISTINCT Polarity FROM Frames")]
        if len(polarities) > 1:
            raise TsfError(
                f"Analysis contains spectra of both polarities {polarities}; "
                "this reader fits one calibration per analysis and does not "
                "support mixed-polarity files.")
        polarity = polarities[0] if polarities else None

        ncal = con.execute(
            "SELECT COUNT(DISTINCT MzCalibration) FROM Frames").fetchone()[0]
        if ncal > 1:
            raise TsfError(
                f"Spectra reference {ncal} different MzCalibration rows; this "
                "reader fits one calibration per analysis and does not support "
                "that.")

        def info(key: str):
            # Prefer the row for this analysis's polarity; accept a row with
            # no polarity recorded; never silently take the other polarity's.
            if polarity is not None:
                row = con.execute(
                    "SELECT Value FROM CalibrationInfo WHERE KeyName=? AND "
                    "(KeyPolarity=? OR KeyPolarity IS NULL) "
                    "ORDER BY KeyPolarity IS NULL",
                    (key, polarity)).fetchone()
            else:
                row = con.execute(
                    "SELECT Value FROM CalibrationInfo WHERE KeyName=?",
                    (key,)).fetchone()
            return row[0] if row else None

        try:
            masses = info("ReferencePeakMasses")
            tofs = info("MeasuredTimesOfFlight")
        except sqlite3.Error as exc:
            raise CalibrationPointsUnavailable(
                f"Cannot read CalibrationInfo: {exc}") from exc

        if masses is None or tofs is None:
            raise CalibrationPointsUnavailable(
                "This analysis records no ReferencePeakMasses / "
                "MeasuredTimesOfFlight in CalibrationInfo, so the calibration "
                "cannot be fitted from the file and m/z cannot be computed.")

        ref_mz = _doubles(masses)
        ref_t = _doubles(tofs)
        n = min(ref_mz.size, ref_t.size)
        ref_mz, ref_t = ref_mz[:n], ref_t[:n]

        ok = np.isfinite(ref_mz) & np.isfinite(ref_t) & (ref_mz > 0) & (ref_t > 0)
        self.ref_mz, self.ref_t = ref_mz[ok], ref_t[ok]
        if self.ref_mz.size < self.order:
            raise CalibrationPointsUnavailable(
                f"Only {self.ref_mz.size} usable calibration points; order "
                f"{self.order} needs at least {self.order}.")

        row = con.execute(
            "SELECT DigitizerTimebase, DigitizerDelay FROM MzCalibration "
            "WHERE Id=(SELECT MzCalibration FROM Frames ORDER BY Id LIMIT 1)"
        ).fetchone()
        if row is None:
            raise CalibrationPointsUnavailable("No MzCalibration row found.")
        self.timebase, self.delay = float(row[0]), float(row[1])
        if self.timebase == 0.0:
            raise TsfError("DigitizerTimebase must be non-zero.")

    # -- the published forward model: m/z -> flight time ------------------

    def _forward(self, p, m):
        m = np.asarray(m, dtype=np.float64)
        t = p[1] + np.sqrt(SCALE / p[0] * np.abs(m))
        if self.order >= 3:
            t = t + p[2] * m
        if self.order >= 4:
            t = t + p[3] * np.abs(m) ** 1.5
        return t

    def _fit(self) -> None:
        m, t = self.ref_mz, self.ref_t

        # Seed from the linear two-point solution: t is affine in sqrt(m).
        A = np.vstack([np.sqrt(m), np.ones_like(m)]).T
        slope, intercept = np.linalg.lstsq(A, t, rcond=None)[0]
        if slope <= 0:
            raise TsfError("Calibration points are not physically ordered.")
        p = np.array([SCALE / slope ** 2, intercept, 0.0, 0.0][:self.order])

        # Gauss-Newton with a numeric Jacobian; 2-4 parameters, converges fast.
        for _ in range(100):
            r = self._forward(p, m) - t
            J = np.empty((m.size, p.size))
            for j in range(p.size):
                h = max(abs(p[j]) * 1e-7, 1e-10)
                pj = p.copy()
                pj[j] += h
                J[:, j] = (self._forward(pj, m) - self._forward(p, m)) / h
            try:
                step, *_ = np.linalg.lstsq(J, -r, rcond=None)
            except np.linalg.LinAlgError:
                break
            p_new = p + step
            if not np.all(np.isfinite(p_new)) or p_new[0] <= 0:
                break
            converged = np.all(np.abs(step) <= 1e-14 * np.maximum(np.abs(p), 1.0))
            p = p_new
            if converged:
                break
        self.params = p

        fitted = self.mz_from_time(t)
        rel = (fitted - m) / m
        self.fit_rms_ppm = float(np.sqrt(np.mean(rel ** 2)) * 1e6)
        self.fit_max_ppm = float(np.max(np.abs(rel)) * 1e6)
        self.mz_min = float(self.ref_mz.min())
        self.mz_max = float(self.ref_mz.max())

    # -- inversion --------------------------------------------------------

    def mz_from_time(self, t) -> np.ndarray:
        """Invert the forward model.

        Flight times at or before the fitted time origin ``ml2`` have no
        physical solution -- no ion arrives before the start of flight -- and
        yield NaN rather than a fabricated value.
        """
        t = np.atleast_1d(np.asarray(t, dtype=np.float64))
        ml1, ml2 = self.params[0], self.params[1]
        d = t - ml2
        valid = d > 0.0
        out = np.full(t.shape, np.nan)
        if not valid.any():
            return out

        # Two-point solution, exact for order 2 and the Newton seed otherwise.
        tv = t[valid]
        m = (d[valid] ** 2) / (SCALE / ml1)
        if self.order == 2:
            out[valid] = m
            return out

        # Newton on the valid entries only, so that convergence is judged on
        # real targets and no divide-by-zero arises from placeholders.
        for _ in range(60):
            f = self._forward(self.params, m) - tv
            h = np.maximum(np.abs(m) * 1e-7, 1e-10)
            fp = (self._forward(self.params, m + h)
                  - self._forward(self.params, m - h)) / (2.0 * h)
            with np.errstate(divide="ignore", invalid="ignore"):
                step = np.where(fp != 0.0, f / fp, 0.0)
            m_new = m - step
            m_new = np.where(m_new > 0.0, m_new, m * 0.5)
            done = np.all(np.abs(m_new - m) <= 1e-13 * np.abs(m))
            m = m_new
            if done:
                break
        out[valid] = m
        return out

    # -- public API -------------------------------------------------------

    def index_to_mz(self, indices):
        scalar = np.isscalar(indices)
        idx = np.atleast_1d(np.asarray(indices, dtype=np.float64))
        t = self.delay + idx * self.timebase
        mz = self.mz_from_time(t)
        if self.warn_outside_range and not self._warned and mz.size:
            finite = np.isfinite(mz)
            if finite.any():
                # Tolerate roundoff at the calibrant endpoints themselves.
                lo = self.mz_min * (1.0 - 1e-9)
                hi = self.mz_max * (1.0 + 1e-9)
                out = (mz[finite] < lo) | (mz[finite] > hi)
                if out.any():
                    self._warned = True   # once per calibration, not per call
                    warnings.warn(
                        f"values outside the calibrant range "
                        f"({self.mz_min:.1f}-{self.mz_max:.1f} m/z) were converted; "
                        "the fitted calibration extrapolates there and accuracy "
                        "degrades (see calibration_fitted docs). Reported once per "
                        "calibration.", RuntimeWarning, stacklevel=2)
        return float(mz[0]) if scalar else mz

    def index_axis(self, n: int) -> np.ndarray:
        """m/z for indices ``0 .. n-1`` -- the axis of a profile spectrum.

        Cached per length: one fit covers the whole analysis, so every
        profile spectrum of the same length shares this axis, and recomputing
        it per spectrum would dominate an imaging run.
        """
        axis = self._axis_cache.get(n)
        if axis is None:
            t = self.delay + np.arange(n, dtype=np.float64) * self.timebase
            axis = self.mz_from_time(t)
            axis.setflags(write=False)
            self._axis_cache[n] = axis
        return axis

    def mz_to_index(self, mzs):
        scalar = np.isscalar(mzs)
        mz = np.atleast_1d(np.asarray(mzs, dtype=np.float64))
        with np.errstate(invalid="ignore"):
            t = np.where(mz >= 0.0, self._forward(self.params, mz), np.nan)
        idx = (t - self.delay) / self.timebase
        return float(idx[0]) if scalar else idx

    def __repr__(self) -> str:
        return (f"FittedCalibration(order={self.order}, "
                f"points={self.ref_mz.size}, "
                f"range={self.mz_min:.1f}-{self.mz_max:.1f}, "
                f"fit_rms={self.fit_rms_ppm:.3f}ppm)")
