#!/usr/bin/env python3
"""Multi-trajectory classical IR spectra from NVE molecular dynamics.

Designed primarily for an isolated water monomer and for i-PI/eT workflows.
The program intentionally applies *no quantum correction factor*.

Main features
-------------
* Pools dipole-derivative correlation functions from independent trajectories
  before Fourier transformation (time-origin or equal-trajectory weighting).
* Accepts scalar polarized dipoles or full Cartesian dipoles.
* Treats i-PI/eT ``L`` entries as photon pseudo-coordinates: they are retained
  in the input trajectory but excluded from every nuclear COM and molecular-
  frame calculation.
* Correctly removes translational dipole current through
      mu_intrinsic = mu_lab - Q * R_COM
  when XYZ coordinates and the molecular net charge are supplied.  For neutral
  H2O, Q=0 and the exact total dipole is translation invariant; COM drift is
  nevertheless diagnosed and the geometry is COM-centered for body-frame work.
* Optional water-fixed frame removes overall molecular rotation and separates
  the in-plane antisymmetric (x), out-of-plane (y), and HOH-bisector (z) axes.
* Five-point dipole differentiation, finite-correlation windowing, explicit
  artificial linewidth, zero padding, trajectory bootstrap confidence bands,
  individual-trajectory spectra, peak tables, and detailed diagnostics.
* No SciPy dependency.  NumPy is required; Matplotlib is optional.

Input dipole text files contain numeric whitespace- or comma-separated columns.
Comments beginning with '#' are ignored.  Automatic column interpretation is:

    4+ columns : time, mux, muy, muz  (first four columns)
    3 columns  : mux, muy, muz        (--dt-fs required)
    2 columns  : time, scalar dipole
    1 column   : scalar dipole        (--dt-fs required)

Override this with --time-col and --dipole-cols (zero-based column indices).
XYZ files are ordinary or extended XYZ trajectories and must contain one frame
per dipole sample.  For more complicated mappings, use a JSON manifest:

[
  {"name": "sim1", "dipole": "Sim1/dipole.dat",
   "xyz": "Sim1/positions.xyz"},
  {"name": "sim2", "dipole": "Sim2/dipole.dat",
   "xyz": "Sim2/positions.xyz"}
]

Paths inside the manifest are resolved relative to the manifest itself.  A
manifest entry may override dt_fs, time_col, dipole_cols, dipole_unit,
time_unit, and net_charge.

Examples
--------
Field-free isotropic spectrum from Cartesian dipoles:

  python water_ir_nve.py --dipoles Sim*/dipole_xyz.dat --dt-fs 0.25 \
      --polarization isotropic --frame lab --output water_free

Laboratory-z cavity spectrum with coordinate diagnostics:

  python water_ir_nve.py --manifest cavity_runs.json --polarization z \
      --frame lab --output cavity_z

Both laboratory and rotation-removed spectra:

  python water_ir_nve.py --manifest runs.json --polarization isotropic \
      --frame both --output water_compare

Run internal numerical tests:

  python water_ir_nve.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import tempfile
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np


VERSION = "1.0.2"
C_CM_PER_FS = 2.99792458e-5
EANGSTROM_TO_AU = 1.8897261254578281
EANGSTROM_TO_DEBYE = 4.803204712570263

ATOMIC_MASSES_U = {
    "H": 1.00784,
    "D": 2.01410177812,
    "T": 3.01604928199,
    "O": 15.999,
}

# i-PI/eT writes the cavity photon coordinate as an extra XYZ pseudo-atom with
# the label "L".  It is not a nucleus, has no atomic mass, and must never enter
# the molecular center of mass or the water-fixed frame.
PHOTON_COORDINATE_LABELS = {"L"}


class InputError(ValueError):
    """Raised for an invalid or inconsistent trajectory input."""


@dataclass
class TrajectorySpec:
    name: str
    dipole: Path
    xyz: Optional[Path] = None
    dt_fs: Optional[float] = None
    time_col: Optional[int] = None
    dipole_cols: Optional[list[int]] = None
    dipole_unit: str = "au"
    time_unit: str = "fs"
    net_charge: float = 0.0


@dataclass
class Diagnostics:
    name: str
    dipole_file: str
    xyz_file: str
    samples_input: int
    samples_used: int
    duration_fs: float
    dt_fs: float
    dt_relative_std: float
    dipole_components: int
    net_charge_e: float
    photon_coordinates_ignored: int
    com_correction_rms_dipole: float
    com_drift_angstrom: float
    com_max_displacement_angstrom: float
    com_mean_speed_angstrom_per_fs: float
    com_max_speed_angstrom_per_fs: float


@dataclass
class PreparedTrajectory:
    spec: TrajectorySpec
    time_fs: np.ndarray
    mu_lab: np.ndarray
    mu_body: Optional[np.ndarray]
    diagnostics: Diagnostics


@dataclass
class CorrelationContribution:
    name: str
    numerator: np.ndarray
    origins: np.ndarray
    acf: np.ndarray


def _positive_float(value: str) -> float:
    x = float(value)
    if not np.isfinite(x) or x <= 0.0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return x


def _fraction(value: str) -> float:
    x = float(value)
    if not 0.0 < x <= 1.0:
        raise argparse.ArgumentTypeError("must lie in (0, 1]")
    return x


def _time_factor_to_fs(unit: str) -> float:
    factors = {"fs": 1.0, "ps": 1000.0, "au": 0.024188843265857}
    try:
        return factors[unit.lower()]
    except KeyError as exc:
        raise InputError(f"unsupported time unit: {unit}") from exc


def _eangstrom_to_dipole_unit(unit: str) -> float:
    factors = {
        "au": EANGSTROM_TO_AU,
        "debye": EANGSTROM_TO_DEBYE,
        "eangstrom": 1.0,
    }
    try:
        return factors[unit.lower()]
    except KeyError as exc:
        raise InputError(f"unsupported dipole unit: {unit}") from exc


def _read_numeric_table(path: Path) -> np.ndarray:
    if not path.is_file():
        raise InputError(f"dipole file not found: {path}")
    rows: list[list[float]] = []
    width: Optional[int] = None
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = line.replace(",", " ").split()
            try:
                row = [float(item) for item in fields]
            except ValueError:
                # Permit one or more un-commented header lines, but never a
                # malformed numeric line after numeric data have started.
                if not rows:
                    continue
                raise InputError(
                    f"non-numeric data in {path} at line {line_number}"
                )
            if width is None:
                width = len(row)
            if len(row) != width:
                raise InputError(
                    f"inconsistent column count in {path} at line {line_number}"
                )
            rows.append(row)
    if not rows:
        raise InputError(f"no numeric data found in {path}")
    data = np.asarray(rows, dtype=float)
    if not np.all(np.isfinite(data)):
        raise InputError(f"NaN or infinite value found in {path}")
    return data


def read_dipole(spec: TrajectorySpec) -> tuple[np.ndarray, np.ndarray, float]:
    data = _read_numeric_table(spec.dipole)
    ncol = data.shape[1]

    if spec.dipole_cols is None:
        if ncol >= 4:
            dipole_cols = [1, 2, 3]
            auto_time_col: Optional[int] = 0
        elif ncol == 3:
            dipole_cols = [0, 1, 2]
            auto_time_col = None
        elif ncol == 2:
            dipole_cols = [1]
            auto_time_col = 0
        else:
            dipole_cols = [0]
            auto_time_col = None
    else:
        dipole_cols = list(spec.dipole_cols)
        auto_time_col = None

    if len(dipole_cols) not in (1, 3):
        raise InputError(
            f"{spec.name}: dipole_cols must select either 1 or 3 columns"
        )
    if min(dipole_cols) < 0 or max(dipole_cols) >= ncol:
        raise InputError(f"{spec.name}: dipole column outside 0..{ncol - 1}")

    time_col = spec.time_col if spec.time_col is not None else auto_time_col
    if time_col is not None:
        if not 0 <= time_col < ncol:
            raise InputError(f"{spec.name}: time column outside 0..{ncol - 1}")
        if time_col in dipole_cols:
            raise InputError(f"{spec.name}: time and dipole columns overlap")
        time_fs = data[:, time_col] * _time_factor_to_fs(spec.time_unit)
        increments = np.diff(time_fs)
        if np.any(increments <= 0.0):
            raise InputError(f"{spec.name}: time values are not strictly increasing")
        dt_fs = float(np.median(increments))
        relative_std = float(np.std(increments) / dt_fs)
        if relative_std > 1.0e-5:
            raise InputError(
                f"{spec.name}: nonuniform time step (relative std={relative_std:.3g}); "
                "resample the dipole series onto a uniform grid first"
            )
        if spec.dt_fs is not None:
            mismatch = abs(dt_fs - spec.dt_fs) / spec.dt_fs
            if mismatch > 1.0e-5:
                raise InputError(
                    f"{spec.name}: --dt-fs={spec.dt_fs:g} conflicts with "
                    f"time-column spacing {dt_fs:g} fs"
                )
    else:
        if spec.dt_fs is None:
            raise InputError(
                f"{spec.name}: no time column detected; supply --dt-fs"
            )
        dt_fs = float(spec.dt_fs)
        relative_std = 0.0
        time_fs = np.arange(data.shape[0], dtype=float) * dt_fs

    if data.shape[0] < 16:
        raise InputError(f"{spec.name}: at least 16 dipole samples are required")
    mu = data[:, dipole_cols].astype(float, copy=True)
    return time_fs, mu, relative_std


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    if not path.is_file():
        raise InputError(f"XYZ file not found: {path}")
    frames: list[np.ndarray] = []
    reference_symbols: Optional[list[str]] = None
    with path.open("r", encoding="utf-8") as handle:
        frame_number = 0
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.strip():
                continue
            frame_number += 1
            try:
                natoms = int(line.strip().split()[0])
            except (ValueError, IndexError) as exc:
                raise InputError(
                    f"invalid atom-count line in {path}, frame {frame_number}"
                ) from exc
            if natoms <= 0:
                raise InputError(f"invalid atom count in {path}, frame {frame_number}")
            if not handle.readline():
                raise InputError(f"missing XYZ comment line in {path}, frame {frame_number}")
            symbols: list[str] = []
            coords = np.empty((natoms, 3), dtype=float)
            for atom in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise InputError(f"truncated {path} at frame {frame_number}")
                fields = atom_line.split()
                if len(fields) < 4:
                    raise InputError(
                        f"invalid XYZ atom line in {path}, frame {frame_number}"
                    )
                symbols.append(fields[0])
                try:
                    coords[atom] = [float(fields[1]), float(fields[2]), float(fields[3])]
                except ValueError as exc:
                    raise InputError(
                        f"invalid XYZ coordinate in {path}, frame {frame_number}"
                    ) from exc
            if reference_symbols is None:
                reference_symbols = symbols
            elif symbols != reference_symbols:
                raise InputError(f"atom ordering changes in {path}, frame {frame_number}")
            frames.append(coords)
    if not frames or reference_symbols is None:
        raise InputError(f"no XYZ frames found in {path}")
    xyz = np.stack(frames)
    if not np.all(np.isfinite(xyz)):
        raise InputError(f"NaN or infinite coordinate found in {path}")
    return reference_symbols, xyz


def atomic_masses(symbols: Sequence[str]) -> np.ndarray:
    masses = []
    for symbol in symbols:
        normalized = symbol.strip().capitalize()
        if normalized not in ATOMIC_MASSES_U:
            raise InputError(
                f"unknown mass for element '{symbol}'; this analyzer supports H/D/T/O"
            )
        masses.append(ATOMIC_MASSES_U[normalized])
    return np.asarray(masses, dtype=float)


def nuclear_coordinates(
    symbols: Sequence[str], xyz: np.ndarray
) -> tuple[list[str], np.ndarray, int]:
    """Remove cavity photon pseudo-atoms from an XYZ trajectory.

    i-PI/eT represents each photon coordinate by an XYZ entry labelled ``L``.
    Such entries are dynamical coordinates but not material particles, so they
    are deliberately excluded from nuclear masses, COM diagnostics, charged-
    system translational correction, and water-body-frame construction.
    """
    photon_mask = np.asarray(
        [symbol.strip().upper() in PHOTON_COORDINATE_LABELS for symbol in symbols],
        dtype=bool,
    )
    nuclear_mask = ~photon_mask
    if not np.any(nuclear_mask):
        raise InputError("XYZ trajectory contains photon coordinates but no nuclei")
    nuclear_symbols = [
        symbol for symbol, keep in zip(symbols, nuclear_mask) if bool(keep)
    ]
    return nuclear_symbols, xyz[:, nuclear_mask, :], int(np.sum(photon_mask))


def center_of_mass(xyz: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.einsum("fai,a->fi", xyz, masses) / float(np.sum(masses))


def water_body_components(
    symbols: Sequence[str], xyz_centered: np.ndarray, mu_lab: np.ndarray
) -> np.ndarray:
    if mu_lab.shape[1] != 3:
        raise InputError("body-frame analysis requires three Cartesian dipole components")
    oxygen = [i for i, s in enumerate(symbols) if s.capitalize() == "O"]
    hydrogens = [i for i, s in enumerate(symbols) if s.capitalize() in {"H", "D", "T"}]
    if len(oxygen) != 1 or len(hydrogens) != 2 or len(symbols) != 3:
        raise InputError(
            "body-frame mode currently requires exactly one O and two H/D/T atoms"
        )
    o = xyz_centered[:, oxygen[0], :]
    r1 = xyz_centered[:, hydrogens[0], :] - o
    r2 = xyz_centered[:, hydrogens[1], :] - o
    norm1 = np.linalg.norm(r1, axis=1)
    norm2 = np.linalg.norm(r2, axis=1)
    if np.any(norm1 < 0.2) or np.any(norm2 < 0.2):
        raise InputError("unphysical O-H distance encountered while building body frame")
    u1 = r1 / norm1[:, None]
    u2 = r2 / norm2[:, None]

    # z: HOH bisector; provisional x: antisymmetric in-plane direction.
    z_axis = u1 + u2
    x_axis = u1 - u2
    z_norm = np.linalg.norm(z_axis, axis=1)
    x_norm = np.linalg.norm(x_axis, axis=1)
    if np.any(z_norm < 1.0e-8) or np.any(x_norm < 1.0e-8):
        raise InputError("degenerate water geometry encountered while building body frame")
    z_axis /= z_norm[:, None]
    x_axis /= x_norm[:, None]
    y_axis = np.cross(z_axis, x_axis)
    y_axis /= np.linalg.norm(y_axis, axis=1)[:, None]
    # Re-orthogonalize x so numerical noise cannot accumulate.
    x_axis = np.cross(y_axis, z_axis)
    x_axis /= np.linalg.norm(x_axis, axis=1)[:, None]
    axes = np.stack((x_axis, y_axis, z_axis), axis=1)
    return np.einsum("fai,fi->fa", axes, mu_lab)


def remove_linear_trend(values: np.ndarray, time_fs: np.ndarray) -> np.ndarray:
    centered_time = time_fs - float(np.mean(time_fs))
    denominator = float(np.dot(centered_time, centered_time))
    if denominator == 0.0:
        return values - np.mean(values, axis=0, keepdims=True)
    slopes = centered_time @ values / denominator
    intercepts = np.mean(values, axis=0)
    return values - intercepts - centered_time[:, None] * slopes


def prepare_trajectory(spec: TrajectorySpec, linear_detrend: bool) -> PreparedTrajectory:
    time_fs, mu, dt_relative_std = read_dipole(spec)
    xyz_file = ""
    mu_body: Optional[np.ndarray] = None
    com_correction_rms = 0.0
    com_drift = float("nan")
    com_max_displacement = float("nan")
    com_mean_speed = float("nan")
    com_max_speed = float("nan")
    photon_coordinates_ignored = 0

    if spec.xyz is not None:
        symbols_all, xyz_all = read_xyz(spec.xyz)
        xyz_file = str(spec.xyz)
        if xyz_all.shape[0] != mu.shape[0]:
            raise InputError(
                f"{spec.name}: {xyz_all.shape[0]} XYZ frames but {mu.shape[0]} dipole samples"
            )
        symbols, xyz, photon_coordinates_ignored = nuclear_coordinates(
            symbols_all, xyz_all
        )
        masses = atomic_masses(symbols)
        r_com = center_of_mass(xyz, masses)
        dt_fs = float(np.median(np.diff(time_fs)))
        v_com = np.gradient(r_com, dt_fs, axis=0, edge_order=2)
        displacement = r_com - r_com[0]
        com_drift = float(np.linalg.norm(displacement[-1]))
        com_max_displacement = float(np.max(np.linalg.norm(displacement, axis=1)))
        speed = np.linalg.norm(v_com, axis=1)
        com_mean_speed = float(np.mean(speed))
        com_max_speed = float(np.max(speed))

        if mu.shape[1] == 3:
            factor = _eangstrom_to_dipole_unit(spec.dipole_unit)
            correction = spec.net_charge * r_com * factor
            mu = mu - correction
            com_correction_rms = float(np.sqrt(np.mean(correction**2)))
            xyz_centered = xyz - r_com[:, None, :]
            mu_body = water_body_components(symbols, xyz_centered, mu)
        elif abs(spec.net_charge) > 1.0e-14:
            raise InputError(
                f"{spec.name}: COM correction for a charged system requires "
                "three Cartesian dipole components"
            )

    if linear_detrend:
        mu = remove_linear_trend(mu, time_fs)
        if mu_body is not None:
            mu_body = remove_linear_trend(mu_body, time_fs)
    else:
        mu = mu - np.mean(mu, axis=0, keepdims=True)
        if mu_body is not None:
            mu_body = mu_body - np.mean(mu_body, axis=0, keepdims=True)

    dt_fs = float(np.median(np.diff(time_fs)))
    diagnostics = Diagnostics(
        name=spec.name,
        dipole_file=str(spec.dipole),
        xyz_file=xyz_file,
        samples_input=int(mu.shape[0]),
        samples_used=max(0, int(mu.shape[0] - 4)),
        duration_fs=float(time_fs[-1] - time_fs[0]),
        dt_fs=dt_fs,
        dt_relative_std=dt_relative_std,
        dipole_components=int(mu.shape[1]),
        net_charge_e=float(spec.net_charge),
        photon_coordinates_ignored=photon_coordinates_ignored,
        com_correction_rms_dipole=com_correction_rms,
        com_drift_angstrom=com_drift,
        com_max_displacement_angstrom=com_max_displacement,
        com_mean_speed_angstrom_per_fs=com_mean_speed,
        com_max_speed_angstrom_per_fs=com_max_speed,
    )
    return PreparedTrajectory(spec, time_fs, mu, mu_body, diagnostics)


def five_point_derivative(values: np.ndarray, dt_fs: float) -> np.ndarray:
    if values.shape[0] < 7:
        raise InputError("five-point differentiation requires at least 7 samples")
    derivative = (
        values[:-4]
        - 8.0 * values[1:-3]
        + 8.0 * values[3:-1]
        - values[4:]
    ) / (12.0 * dt_fs)
    # Remove a residual numerical current; this affects only the zero-frequency bin.
    return derivative - np.mean(derivative, axis=0, keepdims=True)


def select_signal(
    derivative: np.ndarray,
    polarization: str,
    custom_axis: Optional[Sequence[float]],
) -> tuple[np.ndarray, float, str]:
    ncomponents = derivative.shape[1]
    if ncomponents == 1:
        if polarization == "isotropic" or custom_axis is not None:
            raise InputError("isotropic/custom polarization requires Cartesian dipoles")
        return derivative, 1.0, "scalar"

    if custom_axis is not None:
        axis = np.asarray(custom_axis, dtype=float)
        norm = float(np.linalg.norm(axis))
        if norm == 0.0:
            raise InputError("--axis cannot be the zero vector")
        axis /= norm
        projected = derivative @ axis
        label = "axis_" + "_".join(f"{v:.4g}" for v in axis)
        return projected[:, None], 1.0, label

    if polarization == "isotropic":
        # The conventional orientational average is 1/3 of the Cartesian sum.
        return derivative, 1.0 / 3.0, "isotropic"
    component = {"x": 0, "y": 1, "z": 2}[polarization]
    return derivative[:, component : component + 1], 1.0, polarization


def autocorrelation_contribution(
    name: str, signal: np.ndarray, component_factor: float, max_lag: int
) -> CorrelationContribution:
    n = signal.shape[0]
    nfft = 1 << int(math.ceil(math.log2(2 * n - 1)))
    transform = np.fft.rfft(signal, n=nfft, axis=0)
    raw = np.fft.irfft(transform.conj() * transform, n=nfft, axis=0)[:max_lag]
    numerator = component_factor * np.sum(raw, axis=1)
    origins = np.arange(n, n - max_lag, -1, dtype=float)
    acf = numerator / origins
    return CorrelationContribution(name, numerator, origins, acf)


def pool_correlations(
    contributions: Sequence[CorrelationContribution], weighting: str
) -> np.ndarray:
    if not contributions:
        raise InputError("no correlation contributions to pool")
    if weighting == "origins":
        numerator = np.sum([item.numerator for item in contributions], axis=0)
        origins = np.sum([item.origins for item in contributions], axis=0)
        return numerator / origins
    return np.mean([item.acf for item in contributions], axis=0)


def end_taper(n: int, fraction: float) -> np.ndarray:
    taper = np.ones(n, dtype=float)
    if fraction <= 0.0 or n < 2:
        return taper
    start = max(1, int(math.floor((1.0 - fraction) * (n - 1))))
    x = np.linspace(0.0, 1.0, n - start)
    taper[start:] = 0.5 * (1.0 + np.cos(np.pi * x))
    return taper


def correlation_window(
    lag_fs: np.ndarray,
    kind: str,
    broadening_cm1: float,
    taper_fraction: float,
) -> np.ndarray:
    if kind in {"gaussian", "exponential"} and broadening_cm1 <= 0.0:
        raise InputError(f"--broadening-cm1 must be positive for {kind} window")
    if kind == "gaussian":
        # Gaussian spectral FWHM = sqrt(2 ln 2)/(pi*tau*c).
        tau_fs = math.sqrt(2.0 * math.log(2.0)) / (
            math.pi * C_CM_PER_FS * broadening_cm1
        )
        window = np.exp(-0.5 * (lag_fs / tau_fs) ** 2)
    elif kind == "exponential":
        # Lorentzian spectral FWHM = 1/(pi*tau*c).
        tau_fs = 1.0 / (math.pi * C_CM_PER_FS * broadening_cm1)
        window = np.exp(-lag_fs / tau_fs)
    elif kind == "hann":
        x = lag_fs / lag_fs[-1]
        window = 0.5 * (1.0 + np.cos(np.pi * x))
    elif kind == "blackman":
        x = lag_fs / lag_fs[-1]
        window = 0.42 + 0.5 * np.cos(np.pi * x) + 0.08 * np.cos(2.0 * np.pi * x)
    elif kind == "none":
        window = np.ones_like(lag_fs)
    else:
        raise InputError(f"unknown window: {kind}")
    return window * end_taper(len(lag_fs), taper_fraction)


def spectrum_from_acf(
    acf: np.ndarray, dt_fs: float, window: np.ndarray, zero_pad_factor: int
) -> tuple[np.ndarray, np.ndarray]:
    n_lag = len(acf)
    target = max(2, 2 * (n_lag - 1) * zero_pad_factor)
    nfft = 1 << int(math.ceil(math.log2(target)))
    transformed = np.fft.rfft(acf * window, n=nfft)
    # Even extension of the positive-lag correlation: C(0) + 2 sum_{t>0}.
    intensity = dt_fs * (2.0 * transformed.real - acf[0] * window[0])
    frequency_per_fs = np.fft.rfftfreq(nfft, d=dt_fs)
    wavenumber_cm1 = frequency_per_fs / C_CM_PER_FS
    return wavenumber_cm1, intensity


def bootstrap_spectra(
    contributions: Sequence[CorrelationContribution],
    weighting: str,
    dt_fs: float,
    window: np.ndarray,
    zero_pad_factor: int,
    n_bootstrap: int,
    confidence: float,
    seed: int,
    min_cm1: float,
    max_cm1: float,
) -> tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    if n_bootstrap <= 0 or len(contributions) < 2:
        return None, None
    rng = np.random.default_rng(seed)
    spectra = []
    ntraj = len(contributions)
    spectral_mask: Optional[np.ndarray] = None
    for _ in range(n_bootstrap):
        indices = rng.integers(0, ntraj, size=ntraj)
        sampled = [contributions[int(i)] for i in indices]
        acf = pool_correlations(sampled, weighting)
        wavenumber, spectrum = spectrum_from_acf(acf, dt_fs, window, zero_pad_factor)
        if spectral_mask is None:
            spectral_mask = (wavenumber >= min_cm1) & (wavenumber <= max_cm1)
        spectra.append(spectrum[spectral_mask])
    alpha = 0.5 * (100.0 - confidence)
    array = np.asarray(spectra)
    return (
        np.percentile(array, alpha, axis=0),
        np.percentile(array, 100.0 - alpha, axis=0),
    )


def local_peaks(
    x: np.ndarray,
    y: np.ndarray,
    min_x: float,
    max_x: float,
    min_distance: float,
    max_peaks: int,
    min_relative_intensity: float,
) -> list[tuple[float, float]]:
    positive = np.maximum(y, 0.0)
    candidates = np.where(
        (positive[1:-1] > positive[:-2]) & (positive[1:-1] >= positive[2:])
    )[0] + 1
    in_range = (x >= min_x) & (x <= max_x)
    range_maximum = float(np.max(positive[in_range]))
    threshold = min_relative_intensity * range_maximum
    candidates = [
        i for i in candidates
        if min_x <= x[i] <= max_x and positive[i] >= threshold
    ]
    candidates.sort(key=lambda i: positive[i], reverse=True)
    chosen: list[int] = []
    for index in candidates:
        if all(abs(x[index] - x[other]) >= min_distance for other in chosen):
            chosen.append(index)
        if len(chosen) >= max_peaks:
            break
    return sorted([(float(x[i]), float(positive[i])) for i in chosen])


def write_csv(path: Path, header: Sequence[str], columns: Sequence[np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(zip(*columns))


def write_diagnostics(path: Path, diagnostics: Sequence[Diagnostics]) -> None:
    rows = [asdict(item) for item in diagnostics]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: Any) -> Any:
    """Recursively replace non-finite floats so metadata remains strict JSON."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    return value


def plot_spectrum(
    path: Path,
    wavenumber: np.ndarray,
    normalized: np.ndarray,
    ci_low: Optional[np.ndarray],
    ci_high: Optional[np.ndarray],
    min_cm1: float,
    max_cm1: float,
    title: str,
) -> bool:
    try:
        os.environ.setdefault(
            "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "water_ir_matplotlib_cache")
        )
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return False
    mask = (wavenumber >= min_cm1) & (wavenumber <= max_cm1)
    fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
    if ci_low is not None and ci_high is not None:
        ax.fill_between(
            wavenumber[mask], ci_low[mask], ci_high[mask], color="#4C78A8", alpha=0.22,
            linewidth=0.0, label="trajectory bootstrap CI"
        )
    ax.plot(wavenumber[mask], normalized[mask], color="#154c79", linewidth=1.4)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
    ax.set_ylabel("Classical IR intensity / normalized a.u.")
    ax.set_title(title)
    ax.set_xlim(min_cm1, max_cm1)
    ax.set_ylim(bottom=0.0)
    ax.grid(alpha=0.18)
    if ci_low is not None:
        ax.legend(frameon=False)
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return True


def analyze_frame(
    prepared: Sequence[PreparedTrajectory],
    frame: str,
    args: argparse.Namespace,
    output_prefix: Path,
) -> dict[str, Any]:
    derivatives: list[tuple[str, np.ndarray]] = []
    label = ""
    for trajectory in prepared:
        mu = trajectory.mu_lab if frame == "lab" else trajectory.mu_body
        if mu is None:
            raise InputError(
                f"{trajectory.spec.name}: body-frame analysis requested without XYZ coordinates"
            )
        derivative = five_point_derivative(mu, trajectory.diagnostics.dt_fs)
        signal, component_factor, label = select_signal(
            derivative, args.polarization, args.axis
        )
        derivatives.append((trajectory.spec.name, signal))

    shortest = min(signal.shape[0] for _, signal in derivatives)
    if args.max_correlation_fs is not None:
        max_lag = min(shortest, int(math.floor(args.max_correlation_fs / args._dt_fs)) + 1)
    else:
        max_lag = int(math.floor(args.max_correlation_fraction * shortest))
    if max_lag < 8:
        raise InputError("selected correlation interval is too short")

    contributions = [
        autocorrelation_contribution(name, signal, component_factor, max_lag)
        for name, signal in derivatives
    ]
    acf = pool_correlations(contributions, args.weighting)
    lag_fs = np.arange(max_lag, dtype=float) * args._dt_fs
    window = correlation_window(
        lag_fs, args.window, args.broadening_cm1, args.end_taper_fraction
    )
    wavenumber, intensity = spectrum_from_acf(
        acf, args._dt_fs, window, args.zero_pad_factor
    )
    ci_low_raw, ci_high_raw = bootstrap_spectra(
        contributions,
        args.weighting,
        args._dt_fs,
        window,
        args.zero_pad_factor,
        args.bootstrap,
        args.confidence,
        args.seed,
        args.min_cm1,
        args.max_cm1,
    )

    output_mask = (wavenumber >= args.min_cm1) & (wavenumber <= args.max_cm1)
    wavenumber = wavenumber[output_mask]
    intensity = intensity[output_mask]
    normalization_mask = np.isfinite(intensity)
    positive = np.maximum(intensity, 0.0)
    scale = float(np.max(positive[normalization_mask]))
    if scale <= 0.0:
        raise InputError("spectrum has no positive intensity in the requested range")
    normalized = positive / scale
    ci_low = None if ci_low_raw is None else np.maximum(ci_low_raw, 0.0) / scale
    ci_high = None if ci_high_raw is None else np.maximum(ci_high_raw, 0.0) / scale

    tag = f"{frame}_{label}"
    base = output_prefix.parent / f"{output_prefix.name}_{tag}"
    spectrum_columns = [wavenumber, intensity, positive, normalized]
    spectrum_header = [
        "wavenumber_cm-1",
        "intensity_raw",
        "intensity_nonnegative",
        "intensity_normalized",
    ]
    if ci_low is not None and ci_high is not None:
        spectrum_header += [
            f"bootstrap_{args.confidence:g}pct_low_normalized",
            f"bootstrap_{args.confidence:g}pct_high_normalized",
        ]
        spectrum_columns += [ci_low, ci_high]
    write_csv(base.with_name(base.name + "_spectrum.csv"), spectrum_header, spectrum_columns)
    write_csv(
        base.with_name(base.name + "_acf.csv"),
        ["lag_fs", "acf_raw", "acf_normalized", "window", "acf_windowed"],
        [lag_fs, acf, acf / acf[0], window, acf * window],
    )

    individual_files: list[str] = []
    if args.save_individual:
        for contribution in contributions:
            _, individual = spectrum_from_acf(
                contribution.acf, args._dt_fs, window, args.zero_pad_factor
            )
            individual = individual[output_mask]
            safe_name = "".join(
                character if character.isalnum() or character in "-_" else "_"
                for character in contribution.name
            )
            path = base.with_name(base.name + f"_{safe_name}_spectrum.csv")
            write_csv(
                path,
                ["wavenumber_cm-1", "intensity_raw", "intensity_normalized_to_ensemble"],
                [wavenumber, individual, np.maximum(individual, 0.0) / scale],
            )
            individual_files.append(str(path))

    peaks = local_peaks(
        wavenumber,
        normalized,
        args.min_cm1,
        args.max_cm1,
        args.peak_min_distance_cm1,
        args.max_peaks,
        args.peak_min_relative_intensity,
    )
    peak_path = base.with_name(base.name + "_peaks.csv")
    write_csv(
        peak_path,
        ["wavenumber_cm-1", "normalized_intensity"],
        [np.asarray([p[0] for p in peaks]), np.asarray([p[1] for p in peaks])],
    )
    plot_path = base.with_name(base.name + "_spectrum.png")
    plotted = plot_spectrum(
        plot_path,
        wavenumber,
        normalized,
        ci_low,
        ci_high,
        args.min_cm1,
        args.max_cm1,
        f"NVE IR spectrum: {frame} frame, {label}",
    )
    return {
        "frame": frame,
        "polarization": label,
        "correlation_points": max_lag,
        "correlation_time_fs": float(lag_fs[-1]),
        "nominal_resolution_cm-1": float(1.0 / (C_CM_PER_FS * lag_fs[-1])),
        "normalization_scale": scale,
        "peaks": [{"wavenumber_cm-1": p[0], "normalized_intensity": p[1]} for p in peaks],
        "spectrum_file": str(base.with_name(base.name + "_spectrum.csv")),
        "acf_file": str(base.with_name(base.name + "_acf.csv")),
        "peak_file": str(peak_path),
        "plot_file": str(plot_path) if plotted else None,
        "individual_spectrum_files": individual_files,
    }


def load_specs(args: argparse.Namespace) -> list[TrajectorySpec]:
    specs: list[TrajectorySpec] = []
    names_are_explicit = args.manifest is not None
    if args.manifest is not None:
        manifest_path = args.manifest.resolve()
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise InputError(f"manifest not found: {manifest_path}") from exc
        except json.JSONDecodeError as exc:
            raise InputError(f"invalid JSON manifest: {exc}") from exc
        if not isinstance(payload, list) or not payload:
            raise InputError("manifest must be a nonempty JSON list")
        root = manifest_path.parent
        for index, item in enumerate(payload, 1):
            if not isinstance(item, dict) or "dipole" not in item:
                raise InputError(f"manifest entry {index} must contain 'dipole'")
            dipole = Path(item["dipole"])
            xyz = Path(item["xyz"]) if item.get("xyz") else None
            specs.append(
                TrajectorySpec(
                    name=str(item.get("name", f"traj{index}")),
                    dipole=dipole if dipole.is_absolute() else root / dipole,
                    xyz=xyz if xyz is None or xyz.is_absolute() else root / xyz,
                    dt_fs=item.get("dt_fs", args.dt_fs),
                    time_col=item.get("time_col", args.time_col),
                    dipole_cols=item.get("dipole_cols", args.dipole_cols),
                    dipole_unit=str(item.get("dipole_unit", args.dipole_unit)),
                    time_unit=str(item.get("time_unit", args.time_unit)),
                    net_charge=float(item.get("net_charge", args.net_charge)),
                )
            )
    else:
        if not args.dipoles:
            raise InputError("supply --dipoles or --manifest")
        xyz_files = args.xyz or []
        if xyz_files and len(xyz_files) != len(args.dipoles):
            raise InputError("--xyz must contain one file for every --dipoles file")
        # The immediate parent is often identically named across replicas
        # (for example Sim1/Analysis, Sim2/Analysis, ...).  Use the dipole
        # filename as the readable base label and make repeated basenames
        # unique below.  This also prevents --save-individual outputs from
        # overwriting one another.
        used_names: set[str] = set()
        for index, dipole in enumerate(args.dipoles):
            base_name = dipole.stem or f"traj{index + 1}"
            name = base_name
            suffix = 2
            while name in used_names:
                name = f"{base_name}_{suffix}"
                suffix += 1
            used_names.add(name)
            specs.append(
                TrajectorySpec(
                    name=name,
                    dipole=dipole.resolve(),
                    xyz=xyz_files[index].resolve() if xyz_files else None,
                    dt_fs=args.dt_fs,
                    time_col=args.time_col,
                    dipole_cols=args.dipole_cols,
                    dipole_unit=args.dipole_unit,
                    time_unit=args.time_unit,
                    net_charge=args.net_charge,
                )
            )
    names = [spec.name for spec in specs]
    if names_are_explicit and len(set(names)) != len(names):
        raise InputError("trajectory names must be unique")
    return specs


def validate_common_timestep(prepared: Sequence[PreparedTrajectory]) -> float:
    steps = np.asarray([item.diagnostics.dt_fs for item in prepared])
    reference = float(np.median(steps))
    relative = np.max(np.abs(steps - reference)) / reference
    if relative > 1.0e-5:
        raise InputError(
            "trajectories have different time steps; resample them to a common grid "
            f"(found {steps.tolist()} fs)"
        )
    return reference


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Classical multi-trajectory IR spectrum from NVE dipole series",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    source = p.add_argument_group("trajectory input")
    source.add_argument("--dipoles", nargs="+", type=Path, help="dipole text files")
    source.add_argument("--xyz", nargs="+", type=Path, help="matching XYZ trajectories")
    source.add_argument("--manifest", type=Path, help="JSON trajectory manifest")
    source.add_argument("--dt-fs", type=_positive_float, help="sample interval in fs")
    source.add_argument("--time-col", type=int, help="zero-based time column")
    source.add_argument("--dipole-cols", nargs="+", type=int, help="one or three columns")
    source.add_argument("--time-unit", choices=["fs", "ps", "au"], default="fs")
    source.add_argument(
        "--dipole-unit", choices=["au", "debye", "eangstrom"], default="au"
    )
    source.add_argument("--net-charge", type=float, default=0.0, help="molecular charge / e")

    analysis = p.add_argument_group("analysis")
    analysis.add_argument("--frame", choices=["lab", "body", "both"], default="lab")
    analysis.add_argument(
        "--polarization", choices=["x", "y", "z", "isotropic"], default="z"
    )
    analysis.add_argument(
        "--axis", nargs=3, type=float, metavar=("EX", "EY", "EZ"),
        help="custom laboratory/body polarization; overrides --polarization"
    )
    analysis.add_argument(
        "--linear-detrend", action="store_true",
        help="remove linear dipole drift before differentiation"
    )
    analysis.add_argument(
        "--weighting", choices=["origins", "equal"], default="origins",
        help="pool all time origins or give each trajectory equal weight"
    )
    correlation = analysis.add_mutually_exclusive_group()
    correlation.add_argument("--max-correlation-fs", type=_positive_float)
    correlation.add_argument(
        "--max-correlation-fraction", type=_fraction, default=0.5,
        help="fraction of shortest differentiated trajectory"
    )
    analysis.add_argument(
        "--window", choices=["gaussian", "exponential", "hann", "blackman", "none"],
        default="gaussian"
    )
    analysis.add_argument(
        "--broadening-cm1", type=float, default=10.0,
        help="spectral FWHM for Gaussian/exponential windows"
    )
    analysis.add_argument(
        "--end-taper-fraction", type=float, default=0.1,
        help="cosine taper applied to the final fraction of the ACF"
    )
    analysis.add_argument("--zero-pad-factor", type=int, default=4)
    analysis.add_argument("--bootstrap", type=int, default=200)
    analysis.add_argument("--confidence", type=float, default=95.0)
    analysis.add_argument("--seed", type=int, default=1729)

    output = p.add_argument_group("output")
    output.add_argument("--output", type=Path, default=Path("water_ir"))
    output.add_argument("--min-cm1", type=float, default=0.0)
    output.add_argument("--max-cm1", type=float, default=5000.0)
    output.add_argument("--peak-min-distance-cm1", type=float, default=40.0)
    output.add_argument(
        "--peak-min-relative-intensity", type=float, default=0.05,
        help="omit reported peaks below this fraction of the strongest peak"
    )
    output.add_argument("--max-peaks", type=int, default=12)
    output.add_argument("--save-individual", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def validate_args(args: argparse.Namespace) -> None:
    if args.zero_pad_factor < 1:
        raise InputError("--zero-pad-factor must be >= 1")
    if args.bootstrap < 0:
        raise InputError("--bootstrap must be >= 0")
    if not 0.0 < args.confidence < 100.0:
        raise InputError("--confidence must lie in (0, 100)")
    if not 0.0 <= args.end_taper_fraction <= 1.0:
        raise InputError("--end-taper-fraction must lie in [0, 1]")
    if args.min_cm1 < 0.0 or args.max_cm1 <= args.min_cm1:
        raise InputError("invalid --min-cm1/--max-cm1 range")
    if (
        args.max_peaks < 1
        or args.peak_min_distance_cm1 < 0.0
        or not 0.0 <= args.peak_min_relative_intensity <= 1.0
    ):
        raise InputError("invalid peak-selection settings")


def run_self_test() -> None:
    dt = 0.25
    n = 80000  # 20 ps
    time = np.arange(n) * dt
    modes = np.asarray([1595.0, 3657.0, 3756.0])
    amplitudes = np.asarray([1.0, 0.45, 0.40])
    contributions = []
    rng = np.random.default_rng(91)
    max_lag = n // 2
    for trajectory in range(4):
        phases = rng.uniform(0.0, 2.0 * np.pi, size=3)
        mu = sum(
            amplitude * np.cos(2.0 * np.pi * C_CM_PER_FS * mode * time + phase)
            for mode, amplitude, phase in zip(modes, amplitudes, phases)
        )
        mu += rng.normal(0.0, 0.002, size=n)
        derivative = five_point_derivative(mu[:, None], dt)
        contributions.append(
            autocorrelation_contribution(
                f"synthetic{trajectory + 1}", derivative, 1.0, max_lag
            )
        )
    acf = pool_correlations(contributions, "origins")
    lag = np.arange(max_lag) * dt
    window = correlation_window(lag, "gaussian", 8.0, 0.1)
    wn, spectrum = spectrum_from_acf(acf, dt, window, 4)
    for expected in modes:
        neighborhood = (wn > expected - 15.0) & (wn < expected + 15.0)
        observed = float(wn[neighborhood][np.argmax(spectrum[neighborhood])])
        if abs(observed - expected) > 3.0:
            raise AssertionError(f"frequency test failed: {expected} -> {observed}")

    # Charged translating object: mu_lab = mu_intrinsic + Q*R_COM.
    r_com = np.column_stack((0.002 * time, 0.0 * time, 0.0 * time))
    intrinsic = np.column_stack((np.sin(0.01 * time), 0.0 * time, 0.0 * time))
    observed = intrinsic + r_com * EANGSTROM_TO_AU
    recovered = observed - r_com * EANGSTROM_TO_AU
    if not np.allclose(recovered, intrinsic, atol=1.0e-12):
        raise AssertionError("COM correction test failed")

    # Validate frame convention using a rigidly rotated model geometry.
    theta = 0.7
    rotation = np.asarray(
        [[math.cos(theta), -math.sin(theta), 0.0],
         [math.sin(theta), math.cos(theta), 0.0],
         [0.0, 0.0, 1.0]]
    )
    geometry = np.asarray([[0.0, 0.0, 0.0], [0.76, 0.0, 0.59], [-0.76, 0.0, 0.59]])
    geometry = geometry @ rotation.T
    xyz = np.repeat(geometry[None, :, :], 8, axis=0)
    bisector_lab = np.repeat((np.asarray([0.0, 0.0, 1.0]) @ rotation.T)[None, :], 8, axis=0)
    body = water_body_components(["O", "H", "H"], xyz, bisector_lab)
    if not np.allclose(np.abs(body[:, 2]), 1.0, atol=1.0e-12):
        raise AssertionError("body-frame test failed")

    # i-PI/eT cavity trajectories append a massless photon coordinate labelled
    # L.  It must not alter the nuclear COM or prevent body-frame construction.
    xyz_with_photon = np.concatenate(
        (xyz, np.full((xyz.shape[0], 1, 3), 123.456)), axis=1
    )
    nuclear_symbols, nuclear_xyz, photon_count = nuclear_coordinates(
        ["O", "H", "H", "L"], xyz_with_photon
    )
    if photon_count != 1 or nuclear_symbols != ["O", "H", "H"]:
        raise AssertionError("photon-coordinate classification test failed")
    masses = atomic_masses(nuclear_symbols)
    if not np.allclose(
        center_of_mass(nuclear_xyz, masses), center_of_mass(xyz, masses)
    ):
        raise AssertionError("photon coordinate contaminated the nuclear COM")
    photon_body = water_body_components(nuclear_symbols, nuclear_xyz, bisector_lab)
    if not np.allclose(photon_body, body, atol=1.0e-12):
        raise AssertionError("photon coordinate contaminated the water body frame")

    # Exercise manifest parsing, multiple-trajectory pooling, file outputs, and
    # strict-JSON metadata through the same CLI path used in production.
    with tempfile.TemporaryDirectory(prefix="water_ir_selftest_") as temporary:
        root = Path(temporary)

        # Exercise the production input path that combines a scalar polarized
        # dipole with an O/H/H/L i-PI/eT XYZ trajectory.
        cavity_samples = 32
        cavity_dipole = root / "cavity_scalar.dat"
        np.savetxt(cavity_dipole, np.sin(np.arange(cavity_samples) * 0.2))
        cavity_xyz = root / "cavity_positions.xyz"
        xyz_lines: list[str] = []
        for frame_index in range(cavity_samples):
            xyz_lines.extend(
                [
                    "4",
                    f"frame {frame_index}",
                    "O  0.000  0.000  0.000",
                    "H  0.760  0.000  0.590",
                    "H -0.760  0.000  0.590",
                    f"L  0.000  0.000  {10.0 + frame_index:.6f}",
                ]
            )
        cavity_xyz.write_text("\n".join(xyz_lines) + "\n", encoding="utf-8")
        cavity_prepared = prepare_trajectory(
            TrajectorySpec(
                name="cavity_L",
                dipole=cavity_dipole,
                xyz=cavity_xyz,
                dt_fs=dt,
            ),
            linear_detrend=False,
        )
        if cavity_prepared.diagnostics.photon_coordinates_ignored != 1:
            raise AssertionError("O/H/H/L XYZ production-path test failed")

        manifest_entries = []
        n_file = 16000
        for index in range(2):
            phases = rng.uniform(0.0, 2.0 * np.pi, size=3)
            mu_x = sum(
                amplitude * np.cos(
                    2.0 * np.pi * C_CM_PER_FS * mode * time[:n_file] + phase
                )
                for mode, amplitude, phase in zip(modes, amplitudes, phases)
            )
            mu_y = 0.8 * mu_x
            mu_z = 1.2 * mu_x
            table = np.column_stack((time[:n_file], mu_x, mu_y, mu_z))
            dipole_path = root / f"dipole_{index + 1}.dat"
            np.savetxt(dipole_path, table)
            manifest_entries.append(
                {"name": f"test{index + 1}", "dipole": dipole_path.name}
            )
        manifest_path = root / "manifest.json"
        manifest_path.write_text(json.dumps(manifest_entries), encoding="utf-8")
        output = root / "result"
        status = main(
            [
                "--manifest", str(manifest_path),
                "--polarization", "isotropic",
                "--frame", "lab",
                "--bootstrap", "3",
                "--max-cm1", "4500",
                "--output", str(output),
            ]
        )
        if status != 0:
            raise AssertionError("end-to-end CLI test failed")
        metadata_path = root / "result_metadata.json"
        json.loads(metadata_path.read_text(encoding="utf-8"))
        if not (root / "result_lab_isotropic_spectrum.csv").is_file():
            raise AssertionError("end-to-end spectrum output test failed")
    print(
        "Self-test passed: frequencies, COM correction, photon-coordinate "
        "filtering, water body frame, multi-trajectory input, and outputs are valid."
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parser().parse_args(argv)
    if args.self_test:
        run_self_test()
        return 0
    try:
        validate_args(args)
        specs = load_specs(args)
        prepared = [prepare_trajectory(spec, args.linear_detrend) for spec in specs]
        args._dt_fs = validate_common_timestep(prepared)

        if args.frame in {"body", "both"}:
            missing = [item.spec.name for item in prepared if item.mu_body is None]
            if missing:
                raise InputError(
                    "body-frame analysis needs matching XYZ files for: " + ", ".join(missing)
                )
        frames = ["lab", "body"] if args.frame == "both" else [args.frame]
        output_prefix = args.output.resolve()
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        write_diagnostics(
            output_prefix.with_name(output_prefix.name + "_diagnostics.csv"),
            [item.diagnostics for item in prepared],
        )
        results = [analyze_frame(prepared, frame, args, output_prefix) for frame in frames]

        metadata = {
            "program": "water_ir_nve.py",
            "version": VERSION,
            "classical_spectrum": True,
            "quantum_correction_applied": False,
            "dipole_derivative": "five-point centered finite difference",
            "dt_fs": args._dt_fs,
            "trajectory_count": len(prepared),
            "weighting": args.weighting,
            "window": args.window,
            "broadening_cm-1": args.broadening_cm1,
            "end_taper_fraction": args.end_taper_fraction,
            "bootstrap_replicates": args.bootstrap,
            "confidence_percent": args.confidence,
            "com_definition": "mu_intrinsic = mu_lab - net_charge * R_COM",
            "photon_coordinate_note": (
                "XYZ entries labelled L are cavity photon pseudo-coordinates and "
                "are excluded from all nuclear COM and water-body-frame calculations."
            ),
            "neutral_translation_note": (
                "For net_charge=0, the exact total dipole and its derivative are "
                "translation invariant; COM motion is diagnosed but no dipole correction results."
            ),
            "trajectory_diagnostics": [asdict(item.diagnostics) for item in prepared],
            "analyses": results,
        }
        metadata_path = output_prefix.with_name(output_prefix.name + "_metadata.json")
        metadata_path.write_text(
            json.dumps(json_safe(metadata), indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        print(f"Analyzed {len(prepared)} independent trajectory/trajectories at dt={args._dt_fs:g} fs.")
        for result in results:
            print(
                f"  {result['frame']}/{result['polarization']}: correlation time "
                f"{result['correlation_time_fs']:.1f} fs, nominal resolution "
                f"{result['nominal_resolution_cm-1']:.3f} cm^-1"
            )
            for peak in result["peaks"]:
                print(
                    f"    peak {peak['wavenumber_cm-1']:9.3f} cm^-1  "
                    f"I={peak['normalized_intensity']:.4f}"
                )
        print(f"Metadata: {metadata_path}")
        if args.bootstrap > 0 and len(prepared) < 2:
            warnings.warn("bootstrap CI requires at least two trajectories; CI omitted")
        return 0
    except (InputError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
