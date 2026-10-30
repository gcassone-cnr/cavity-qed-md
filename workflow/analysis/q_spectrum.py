#!/usr/bin/env python3
"""Calculate the photon-coordinate spectrum from i-PI XYZ trajectories.

In the eT/i-PI CBO-HF setup used for the water simulations, the auxiliary
particle labelled ``L`` stores the photon displacement in its x coordinate:

    L   q   0.0   0.0

The script extracts q(t), subtracts its mean, evaluates the autocorrelation

    Cqq(t) = <delta q(0) delta q(t)>

over all available time origins, and Fourier transforms the windowed
autocorrelation to obtain

    Sqq(nu) = 2 integral_0^infinity Cqq(t) cos(2 pi nu t) dt.

The spectrum is a photon-coordinate power spectrum.  Unlike the dipole IR
spectrum, it is NOT multiplied by frequency squared.

Examples
--------
Classical, single-bead CBO-HF trajectory:

    python3 q_spectrum.py \
      runs/cavity/h2o_cbohf_ir_nve.positions_au_0.xyz \
      --sample-dt-fs 0.05 \
      --min-wavenumber 3500 --max-wavenumber 4300 \
      --cavity-omega-au 0.0177243 \
      --output-prefix runs/cavity/q_spectrum

If the XYZ trajectory was written every five MD steps and the MD timestep was
0.05 fs, use ``--sample-dt-fs 0.25``.

Several bead files can be supplied with a shell wildcard:

    python3 q_spectrum.py runs/pimd/beads_au_*.xyz \
      --sample-dt-fs 0.1 --output-prefix runs/pimd/q_centroid_spectrum

By default, corresponding q values are averaged over the supplied files before
the autocorrelation is calculated.  This is a bead-centroid diagnostic.  Plain
PIMD trajectories do not, by themselves, provide a physical real-time quantum
spectrum.

Outputs are whitespace-delimited and can be plotted directly with gnuplot:

    plot 'q_spectrum.dat' using 1:5 with lines title 'q spectrum'

Only NumPy is required.  Matplotlib is optional.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Sequence

import numpy as np


C_CM_PER_S = 2.99792458e10
HARTREE_TO_WAVENUMBER = 219474.6313632


@dataclass
class QTrajectory:
    """Photon-coordinate trajectory extracted from one XYZ file."""

    q: np.ndarray
    steps: np.ndarray | None


def _parse_step(comment: str) -> int | None:
    """Extract an i-PI step number from an XYZ comment when present."""
    patterns = (
        r"\bStep\s*[:=]\s*(-?\d+)",
        r"\bstep\s*[:=]\s*(-?\d+)",
        r"\bi\s*[:=]\s*(-?\d+)",
    )
    for pattern in patterns:
        match = re.search(pattern, comment)
        if match:
            return int(match.group(1))
    return None


def read_ipi_xyz(path: Path, label: str, component: str) -> QTrajectory:
    """Read q from the selected labelled particle in every XYZ frame."""
    if not path.is_file():
        raise FileNotFoundError(path)

    component_index = {"x": 1, "y": 2, "z": 3}[component]
    q_values: list[float] = []
    steps: list[int | None] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        frame = 0
        while True:
            first_line = handle.readline()
            if not first_line:
                break
            if not first_line.strip():
                continue
            try:
                natoms = int(first_line.split()[0])
            except (ValueError, IndexError) as exc:
                raise ValueError(
                    f"{path}: expected an XYZ atom count at frame {frame + 1}, "
                    f"but found {first_line.strip()!r}"
                ) from exc
            if natoms <= 0:
                raise ValueError(
                    f"{path}: invalid atom count {natoms} at frame {frame + 1}"
                )

            comment = handle.readline()
            if not comment:
                raise ValueError(
                    f"{path}: missing comment line at frame {frame + 1}"
                )

            selected: list[float] = []
            for atom_index in range(natoms):
                atom_line = handle.readline()
                if not atom_line:
                    raise ValueError(
                        f"{path}: truncated frame {frame + 1}; expected "
                        f"{natoms} atom records"
                    )
                fields = atom_line.split()
                if len(fields) < 4:
                    raise ValueError(
                        f"{path}: malformed atom record at frame {frame + 1}, "
                        f"atom {atom_index + 1}: {atom_line.strip()!r}"
                    )
                if fields[0] == label:
                    try:
                        selected.append(float(fields[component_index]))
                    except ValueError as exc:
                        raise ValueError(
                            f"{path}: non-numerical {component} coordinate for "
                            f"{label!r} at frame {frame + 1}"
                        ) from exc

            if len(selected) != 1:
                raise ValueError(
                    f"{path}: expected exactly one atom labelled {label!r} in "
                    f"frame {frame + 1}, found {len(selected)}"
                )
            q_values.append(selected[0])
            steps.append(_parse_step(comment))
            frame += 1

    if len(q_values) < 4:
        raise ValueError(
            f"{path}: only {len(q_values)} frames were found; at least four "
            "are required"
        )
    q = np.asarray(q_values, dtype=float)
    if not np.all(np.isfinite(q)):
        raise ValueError(f"{path}: non-finite photon coordinate found")

    if all(step is not None for step in steps):
        step_array: np.ndarray | None = np.asarray(steps, dtype=int)
        if np.any(np.diff(step_array) <= 0):
            raise ValueError(f"{path}: XYZ step numbers are not increasing")
    else:
        step_array = None
    return QTrajectory(q=q, steps=step_array)


def load_file_average(
    paths: Sequence[Path], label: str, component: str
) -> tuple[np.ndarray, np.ndarray | None]:
    """Average corresponding q trajectories over files/beads."""
    if not paths:
        raise ValueError("At least one XYZ trajectory is required")

    trajectories = [read_ipi_xyz(path, label, component) for path in paths]
    lengths = {len(trajectory.q) for trajectory in trajectories}
    if len(lengths) != 1:
        raise ValueError(
            "The supplied XYZ files contain different numbers of frames: "
            f"{sorted(lengths)}"
        )

    reference_steps = trajectories[0].steps
    for trajectory in trajectories[1:]:
        if (reference_steps is None) != (trajectory.steps is None):
            raise ValueError(
                "Step numbers could be parsed from only some of the XYZ files"
            )
        if reference_steps is not None and not np.array_equal(
            reference_steps, trajectory.steps
        ):
            raise ValueError(
                "The supplied XYZ files do not contain matching step arrays"
            )

    stacked = np.stack([trajectory.q for trajectory in trajectories], axis=0)
    return np.mean(stacked, axis=0), reference_steps


def autocorrelation_all_origins(signal: np.ndarray) -> np.ndarray:
    """Return the unbiased linear autocorrelation using an FFT."""
    nframes = len(signal)
    nfft = 1 << (2 * nframes - 1).bit_length()
    transformed = np.fft.rfft(signal, n=nfft)
    correlation = np.fft.irfft(
        transformed.conj() * transformed, n=nfft
    )[:nframes]
    correlation /= np.arange(nframes, 0, -1)
    return correlation


def correlation_window(length: int, mode: str) -> np.ndarray:
    """Return a one-sided window that is one at zero lag."""
    if length < 2 or mode == "none":
        return np.ones(length)
    if mode == "hann":
        return np.hanning(2 * length - 1)[length - 1 :]
    if mode == "blackman":
        return np.blackman(2 * length - 1)[length - 1 :]
    raise ValueError(f"Unknown correlation window: {mode}")


def calculate_spectrum(
    q: np.ndarray,
    sample_dt_fs: float,
    discard_frames: int,
    max_correlation_frames: int | None,
    window_mode: str,
    zero_padding: int,
) -> dict[str, np.ndarray]:
    """Calculate q fluctuations, Cqq(t), and Sqq(nu)."""
    if sample_dt_fs <= 0.0:
        raise ValueError("--sample-dt-fs must be positive")
    if zero_padding < 1:
        raise ValueError("--zero-padding must be at least one")
    if not 0 <= discard_frames < len(q) - 3:
        raise ValueError("Discard interval leaves fewer than four samples")

    retained_q = np.asarray(q[discard_frames:], dtype=float)
    q_mean = float(np.mean(retained_q))
    fluctuations = retained_q - q_mean
    acf = autocorrelation_all_origins(fluctuations)

    if max_correlation_frames is None:
        correlation_length = len(acf)
    else:
        correlation_length = min(max_correlation_frames, len(acf))
    if correlation_length < 4:
        raise ValueError("Correlation interval contains fewer than four samples")

    acf = acf[:correlation_length]
    time_fs = np.arange(correlation_length) * sample_dt_fs
    window = correlation_window(correlation_length, window_mode)
    windowed_acf = acf * window

    transform_length = zero_padding * correlation_length
    fft_acf = np.real(np.fft.rfft(windowed_acf, n=transform_length))
    # Correct one-sided cosine transform: the zero-lag term occurs only once.
    spectral_density = sample_dt_fs * (
        2.0 * fft_acf - windowed_acf[0]
    )

    frequency_hz = np.fft.rfftfreq(
        transform_length, d=sample_dt_fs * 1.0e-15
    )
    wavenumber = frequency_hz / C_CM_PER_S
    frequency_thz = frequency_hz * 1.0e-12

    return {
        "retained_q": retained_q,
        "fluctuations": fluctuations,
        "q_mean": np.asarray([q_mean]),
        "time_series_fs": np.arange(len(retained_q)) * sample_dt_fs,
        "acf_time_fs": time_fs,
        "acf": acf,
        "window": window,
        "windowed_acf": windowed_acf,
        "wavenumber_cm-1": wavenumber,
        "frequency_thz": frequency_thz,
        "spectral_density": spectral_density,
    }


def normalize_spectrum(
    wavenumber: np.ndarray,
    spectral_density: np.ndarray,
    minimum: float,
    maximum: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Clip finite-sampling negatives and normalize in a selected interval."""
    if maximum <= minimum:
        raise ValueError("--max-wavenumber must exceed --min-wavenumber")
    selected = (wavenumber >= minimum) & (wavenumber <= maximum)
    if not np.any(selected):
        raise ValueError("Requested wavenumber interval contains no FFT points")

    nonnegative = np.maximum(spectral_density, 0.0)
    scale = float(np.max(nonnegative[selected]))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(
            "The q spectrum has no positive intensity in the requested interval"
        )
    return nonnegative, nonnegative / scale


def find_peaks(
    x: np.ndarray,
    y: np.ndarray,
    minimum: float,
    maximum: float,
    threshold: float,
    minimum_separation: float,
    maximum_peaks: int,
) -> list[tuple[float, float]]:
    """Find the strongest separated local maxima without requiring SciPy."""
    local = np.flatnonzero(
        (y[1:-1] > y[:-2])
        & (y[1:-1] >= y[2:])
        & (y[1:-1] >= threshold)
    ) + 1
    local = local[(x[local] >= minimum) & (x[local] <= maximum)]
    ranked = local[np.argsort(y[local])[::-1]]

    chosen: list[int] = []
    for index in ranked:
        if all(
            abs(x[index] - x[other]) >= minimum_separation
            for other in chosen
        ):
            chosen.append(int(index))
        if len(chosen) == maximum_peaks:
            break
    return sorted(
        ((float(x[index]), float(y[index])) for index in chosen),
        key=lambda pair: pair[0],
    )


def save_outputs(
    output_prefix: Path,
    result: dict[str, np.ndarray],
    nonnegative: np.ndarray,
    normalized: np.ndarray,
    minimum: float,
    maximum: float,
    cavity_wavenumber: float | None,
) -> tuple[Path, Path, Path, Path | None]:
    """Write whitespace-delimited data files and an optional PNG plot."""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    spectrum_path = output_prefix.with_suffix(".dat")
    acf_path = output_prefix.with_name(output_prefix.name + "_acf").with_suffix(
        ".dat"
    )
    q_path = output_prefix.with_name(output_prefix.name + "_q").with_suffix(
        ".dat"
    )
    plot_path = output_prefix.with_suffix(".png")

    np.savetxt(
        spectrum_path,
        np.column_stack(
            (
                result["wavenumber_cm-1"],
                result["frequency_thz"],
                result["spectral_density"],
                nonnegative,
                normalized,
            )
        ),
        fmt="%.12e",
        header=(
            "Columns:\n"
            "1 wavenumber_cm-1\n"
            "2 frequency_THz\n"
            "3 Sqq_raw_au2_fs\n"
            "4 Sqq_nonnegative_au2_fs\n"
            "5 Sqq_normalized"
        ),
        comments="# ",
    )

    acf0 = float(result["acf"][0])
    normalized_acf = (
        result["acf"] / acf0 if abs(acf0) > 0.0 else np.zeros_like(result["acf"])
    )
    np.savetxt(
        acf_path,
        np.column_stack(
            (
                result["acf_time_fs"],
                result["acf"],
                normalized_acf,
                result["window"],
                result["windowed_acf"],
            )
        ),
        fmt="%.12e",
        header=(
            "Columns:\n"
            "1 correlation_time_fs\n"
            "2 Cqq_au2\n"
            "3 Cqq_normalized\n"
            "4 window\n"
            "5 windowed_Cqq_au2"
        ),
        comments="# ",
    )

    np.savetxt(
        q_path,
        np.column_stack(
            (
                result["time_series_fs"],
                result["retained_q"],
                result["fluctuations"],
            )
        ),
        fmt="%.12e",
        header=(
            "Columns:\n"
            "1 time_fs\n"
            "2 q_au\n"
            "3 q_minus_mean_au"
        ),
        comments="# ",
    )

    written_plot: Path | None = plot_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "WARNING: matplotlib is unavailable; data files were written but "
            "the PNG plot was skipped.",
            file=sys.stderr,
        )
        written_plot = None
    else:
        selected = (
            (result["wavenumber_cm-1"] >= minimum)
            & (result["wavenumber_cm-1"] <= maximum)
        )
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(
            result["wavenumber_cm-1"][selected],
            normalized[selected],
            color="#6a3d9a",
            linewidth=1.5,
            label=r"$S_{qq}$",
        )
        if cavity_wavenumber is not None:
            ax.axvline(
                cavity_wavenumber,
                color="#d95f02",
                linestyle="--",
                linewidth=1.1,
                label=(
                    rf"$\omega_c={cavity_wavenumber:.1f}\ "
                    r"\mathrm{cm}^{-1}$"
                ),
            )
        ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
        ax.set_ylabel(r"Normalized photon spectrum $S_{qq}$")
        ax.set_xlim(minimum, maximum)
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.2)
        if cavity_wavenumber is not None:
            ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=220)
        plt.close(fig)

    return spectrum_path, acf_path, q_path, written_plot


def run_self_test() -> int:
    """Verify recovery of two synthetic photon-coordinate frequencies."""
    targets = np.asarray([3850.0, 3965.0])
    sample_dt_fs = 0.25
    time_fs = np.arange(0.0, 15000.0, sample_dt_fs)
    q = np.zeros_like(time_fs)
    for amplitude, target in zip((1.0, 0.7), targets):
        cycles_per_fs = target * C_CM_PER_S * 1.0e-15
        q += amplitude * np.sin(2.0 * np.pi * cycles_per_fs * time_fs)

    result = calculate_spectrum(
        q=q,
        sample_dt_fs=sample_dt_fs,
        discard_frames=0,
        max_correlation_frames=None,
        window_mode="hann",
        zero_padding=4,
    )
    _, normalized = normalize_spectrum(
        result["wavenumber_cm-1"],
        result["spectral_density"],
        3700.0,
        4100.0,
    )
    peaks = find_peaks(
        result["wavenumber_cm-1"],
        normalized,
        3700.0,
        4100.0,
        threshold=0.1,
        minimum_separation=40.0,
        maximum_peaks=2,
    )
    recovered = np.asarray([peak[0] for peak in peaks])
    if len(recovered) != 2 or np.max(np.abs(recovered - targets)) > 2.0:
        raise AssertionError(
            f"Self-test recovered {recovered.tolist()} cm^-1; expected "
            f"{targets.tolist()} cm^-1"
        )
    print(
        "Self-test passed: recovered peaks at "
        + ", ".join(f"{value:.3f} cm^-1" for value in recovered)
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the photon coordinate from i-PI XYZ trajectories and "
            "calculate its autocorrelation spectrum."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "xyz_files",
        nargs="*",
        type=Path,
        help=(
            "i-PI XYZ trajectory file(s). Multiple files are interpreted as "
            "corresponding bead trajectories and averaged before analysis."
        ),
    )
    parser.add_argument(
        "--sample-dt-fs",
        type=float,
        help=(
            "Time in fs between saved XYZ frames: MD timestep times XYZ "
            "trajectory stride. Required unless --self-test is used."
        ),
    )
    parser.add_argument(
        "--label",
        default="L",
        help="XYZ particle label whose coordinate stores q.",
    )
    parser.add_argument(
        "--component",
        choices=["x", "y", "z"],
        default="x",
        help="Coordinate of the selected particle that stores q.",
    )
    parser.add_argument(
        "--discard-fs",
        type=float,
        default=0.0,
        help="Initial retained-trajectory interval to discard.",
    )
    parser.add_argument(
        "--max-correlation-fs",
        type=float,
        help=(
            "Maximum correlation time. The complete retained trajectory is "
            "used when omitted."
        ),
    )
    parser.add_argument(
        "--window",
        choices=["hann", "blackman", "none"],
        default="hann",
        help="Window applied to the one-sided autocorrelation.",
    )
    parser.add_argument(
        "--zero-padding",
        type=int,
        default=4,
        help="FFT zero-padding factor; interpolation only, not resolution.",
    )
    parser.add_argument(
        "--min-wavenumber",
        type=float,
        default=0.0,
        help="Lower limit for normalization, plotting and peak finding.",
    )
    parser.add_argument(
        "--max-wavenumber",
        type=float,
        default=5000.0,
        help="Upper limit for normalization, plotting and peak finding.",
    )
    parser.add_argument(
        "--cavity-omega-au",
        type=float,
        help=(
            "Optional cavity angular frequency in atomic units. Its "
            "wavenumber is reported and marked in the plot."
        ),
    )
    parser.add_argument(
        "--peak-threshold",
        type=float,
        default=0.05,
        help="Minimum normalized Sqq reported as a peak.",
    )
    parser.add_argument(
        "--peak-separation",
        type=float,
        default=25.0,
        help="Minimum separation between reported peaks in cm^-1.",
    )
    parser.add_argument(
        "--max-peaks",
        type=int,
        default=10,
        help="Maximum number of peaks printed.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("q_spectrum"),
        help="Output path prefix.",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run a synthetic two-frequency test without reading an XYZ file.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.self_test:
        return run_self_test()
    if not args.xyz_files:
        raise ValueError("Provide at least one i-PI XYZ trajectory")
    if args.sample_dt_fs is None or args.sample_dt_fs <= 0.0:
        raise ValueError(
            "--sample-dt-fs is required and must be positive. It is the time "
            "between two saved XYZ frames."
        )
    if args.discard_fs < 0.0:
        raise ValueError("--discard-fs cannot be negative")
    if args.max_correlation_fs is not None and args.max_correlation_fs <= 0.0:
        raise ValueError("--max-correlation-fs must be positive")
    if args.cavity_omega_au is not None and args.cavity_omega_au <= 0.0:
        raise ValueError("--cavity-omega-au must be positive")

    q, steps = load_file_average(
        args.xyz_files, args.label, args.component
    )
    discard_frames = int(round(args.discard_fs / args.sample_dt_fs))
    maximum_frames = (
        None
        if args.max_correlation_fs is None
        else int(round(args.max_correlation_fs / args.sample_dt_fs)) + 1
    )
    result = calculate_spectrum(
        q=q,
        sample_dt_fs=args.sample_dt_fs,
        discard_frames=discard_frames,
        max_correlation_frames=maximum_frames,
        window_mode=args.window,
        zero_padding=args.zero_padding,
    )
    nonnegative, normalized = normalize_spectrum(
        result["wavenumber_cm-1"],
        result["spectral_density"],
        args.min_wavenumber,
        args.max_wavenumber,
    )
    cavity_wavenumber = (
        None
        if args.cavity_omega_au is None
        else args.cavity_omega_au * HARTREE_TO_WAVENUMBER
    )
    spectrum_path, acf_path, q_path, plot_path = save_outputs(
        args.output_prefix,
        result,
        nonnegative,
        normalized,
        args.min_wavenumber,
        args.max_wavenumber,
        cavity_wavenumber,
    )
    peaks = find_peaks(
        result["wavenumber_cm-1"],
        normalized,
        args.min_wavenumber,
        args.max_wavenumber,
        args.peak_threshold,
        args.peak_separation,
        args.max_peaks,
    )

    saved_stride = "unavailable"
    if steps is not None and len(steps) > 1:
        differences = np.diff(steps)
        if np.all(differences == differences[0]):
            saved_stride = str(int(differences[0]))

    retained_frames = len(q) - discard_frames
    duration_fs = max(retained_frames - 1, 0) * args.sample_dt_fs
    correlation_duration_fs = float(result["acf_time_fs"][-1])
    nominal_resolution = (
        1.0 / (C_CM_PER_S * correlation_duration_fs * 1.0e-15)
        if correlation_duration_fs > 0.0
        else float("inf")
    )

    print(f"XYZ files averaged: {len(args.xyz_files)}")
    print(f"Particle/component: {args.label}.{args.component}")
    print(f"Frames retained: {retained_frames}")
    print(f"Saved-step stride reported by i-PI: {saved_stride}")
    print(f"Sample interval: {args.sample_dt_fs:.8g} fs")
    print(f"Retained duration: {duration_fs:.3f} fs")
    print(f"Mean q: {float(result['q_mean'][0]):.10e} au")
    print(f"Standard deviation q: {np.std(result['retained_q']):.10e} au")
    print(f"Correlation duration: {correlation_duration_fs:.3f} fs")
    print(f"Nominal Fourier resolution: {nominal_resolution:.3f} cm^-1")
    if cavity_wavenumber is not None:
        print(
            f"Bare cavity frequency: {cavity_wavenumber:.3f} cm^-1 "
            f"({args.cavity_omega_au:.10g} au)"
        )
    print(f"Spectrum: {spectrum_path}")
    print(f"Autocorrelation: {acf_path}")
    print(f"Extracted q trajectory: {q_path}")
    if plot_path is not None:
        print(f"Plot: {plot_path}")
    if peaks:
        print("Strongest separated q-spectrum peaks:")
        for wavenumber, intensity in peaks:
            print(f"  {wavenumber:10.3f} cm^-1  {intensity:8.5f}")
    else:
        print("No local maxima passed the requested peak threshold.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
