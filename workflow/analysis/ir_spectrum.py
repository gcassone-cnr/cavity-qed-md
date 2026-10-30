#!/usr/bin/env python3
"""Calculate an IR spectrum from i-PI dipole ``extras`` trajectories.

The script reads the dipole produced by an i-PI trajectory such as

    <trajectory filename='dipole_z_au'
                stride='5'
                extra_type='polarized_dipole_au'>
      extras
    </trajectory>

For a classical one-bead trajectory, pass the single ``dipole_z_au`` file.
If several bead files are supplied, the dipole is averaged over the beads
before the correlation function is evaluated, as required for RPMD/TRPMD.

The reported relative IR intensity is

    I(w) proportional to w**2 FourierTransform[<dmu(0)dmu(t)>],

where dmu is the fluctuation of the selected dipole component.  The overall
temperature-dependent prefactor is irrelevant after normalization.

Example
-------
Run with a Python environment containing NumPy (and optionally Matplotlib):

    python3 ir_spectrum.py \
      runs/ir_nve_z/h2o_cbohf_ir_nve.dipole_z_au_0 \
      --sample-dt-fs 0.25 \
      --output-prefix runs/ir_nve_z/ir_z

This writes ``ir_z.csv``, ``ir_z_acf.csv`` and ``ir_z.png``.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Sequence

import numpy as np


C_CM_PER_S = 2.99792458e10


def read_ipi_extras(
    path: Path, extra_key: str | None
) -> tuple[np.ndarray, np.ndarray | None, str]:
    """Read a numerical i-PI extras file without requiring ASE.

    i-PI's public ``read_trajectory`` helper currently checks for ASE before
    dispatching to the extras parser.  The extras format itself is simple and
    does not contain atomic structures, so it is parsed directly here:

        #EXTRAS(property_name)# Step: step_number
        numerical data
    """
    header_pattern = re.compile(
        r"^\s*#EXTRAS[^(]*\(([^)]+)\)#\s*Step:\s*(\d+)"
    )
    property_name: str | None = None
    steps: list[int] = []
    frames: list[np.ndarray] = []
    frame_lines: list[str] = []

    def finish_frame() -> None:
        if not frame_lines:
            return
        text = " ".join(line.strip() for line in frame_lines if line.strip())
        values = np.fromstring(text, sep=" ", dtype=float)
        if values.size == 0:
            raise ValueError(
                f"Could not parse numerical extras data near step "
                f"{steps[-1] if steps else 'unknown'} in {path}"
            )
        frames.append(values)
        frame_lines.clear()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = header_pattern.match(line)
            if match:
                finish_frame()
                current_property = match.group(1).strip()
                if property_name is None:
                    property_name = current_property
                elif current_property != property_name:
                    raise ValueError(
                        f"Inconsistent extras properties in {path}: "
                        f"{property_name!r} and {current_property!r}"
                    )
                steps.append(int(match.group(2)))
            elif property_name is not None and line.strip():
                frame_lines.append(line)
    finish_frame()

    if property_name is None or not steps:
        raise ValueError(
            f"No '#EXTRAS(property)# Step:' records were found in {path}"
        )
    if len(frames) != len(steps):
        raise ValueError(
            f"{path} contains {len(steps)} step headers but "
            f"{len(frames)} numerical frames"
        )
    frame_sizes = {frame.size for frame in frames}
    if len(frame_sizes) != 1:
        raise ValueError(
            f"Inconsistent numerical frame sizes in {path}: {frame_sizes}"
        )
    if extra_key is not None and extra_key != property_name:
        raise KeyError(
            f"Requested extras field {extra_key!r}, but {path} contains "
            f"{property_name!r}"
        )

    values = np.asarray(frames, dtype=float)
    values = np.squeeze(values)
    if values.ndim == 0:
        values = values.reshape(1)
    if values.ndim not in (1, 2):
        raise ValueError(
            f"Dipole data in {path} have unsupported shape {values.shape}"
        )
    if values.ndim == 2 and values.shape[1] == 1:
        values = values[:, 0]
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Non-finite dipole value found in {path}")

    step_array = np.asarray(steps, dtype=float)
    if len(step_array) != len(values):
        raise ValueError(
            f"Step and dipole arrays have different lengths in {path}"
        )
    return values, step_array, property_name


def select_components(values: np.ndarray, component: str) -> np.ndarray:
    """Return an (nframes, nchannels) array for the requested observable."""
    if values.ndim == 1:
        if component not in {"auto", "polarized"}:
            raise ValueError(
                "A scalar polarized-dipole file was supplied, so --component "
                "must be 'auto' or 'polarized'."
            )
        return values[:, None]

    if values.shape[1] != 3:
        raise ValueError(
            f"Vector dipoles must have three columns; received {values.shape}"
        )
    if component == "auto":
        component = "isotropic"
    if component == "isotropic":
        # The isotropic spectrum is the mean of the x, y and z correlations.
        return values
    if component in {"x", "y", "z"}:
        return values[:, "xyz".index(component)][:, None]
    raise ValueError(
        "--component polarized requires a scalar polarized-dipole trajectory"
    )


def load_bead_average(
    paths: Sequence[Path], extra_key: str | None, component: str
) -> tuple[np.ndarray, np.ndarray | None, str]:
    """Read corresponding bead files and form the bead-averaged dipole."""
    if not paths:
        raise ValueError("At least one dipole file is required")

    arrays: list[np.ndarray] = []
    reference_steps: np.ndarray | None = None
    selected_keys: list[str] = []
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        values, steps, selected_key = read_ipi_extras(path, extra_key)
        values = select_components(values, component)
        arrays.append(values)
        selected_keys.append(selected_key)
        if reference_steps is None and steps is not None:
            reference_steps = steps
        elif steps is not None and reference_steps is not None:
            if len(steps) != len(reference_steps) or not np.allclose(
                steps, reference_steps, rtol=0.0, atol=1.0e-12
            ):
                raise ValueError(
                    "The supplied bead files do not contain matching step arrays"
                )

    shapes = {array.shape for array in arrays}
    if len(shapes) != 1:
        raise ValueError(f"Bead dipole arrays have inconsistent shapes: {shapes}")
    if len(set(selected_keys)) != 1:
        raise ValueError(
            f"Different extras fields were selected: {sorted(set(selected_keys))}"
        )
    return np.mean(np.stack(arrays, axis=0), axis=0), reference_steps, selected_keys[0]


def autocorrelation_multiple_origins(signal: np.ndarray) -> np.ndarray:
    """Unbiased linear autocorrelation, averaged over all time origins."""
    nframes = signal.shape[0]
    nfft = 1 << (2 * nframes - 1).bit_length()
    transformed = np.fft.rfft(signal, n=nfft, axis=0)
    correlation = np.fft.irfft(
        transformed.conj() * transformed, n=nfft, axis=0
    )[:nframes]
    correlation /= np.arange(nframes, 0, -1)[:, None]
    return np.mean(correlation, axis=1)


def correlation_window(length: int, mode: str) -> np.ndarray:
    """Return a one-sided window equal to one at zero lag."""
    if length < 2 or mode == "none":
        return np.ones(length)
    if mode == "hann":
        return np.hanning(2 * length - 1)[length - 1 :]
    if mode == "blackman":
        return np.blackman(2 * length - 1)[length - 1 :]
    raise ValueError(f"Unknown window: {mode}")


def calculate_ir(
    dipole: np.ndarray,
    sample_dt_fs: float,
    discard_frames: int,
    max_correlation_frames: int | None,
    window_mode: str,
    zero_padding: int,
) -> dict[str, np.ndarray]:
    """Calculate the autocorrelation and relative IR spectrum."""
    if sample_dt_fs <= 0.0:
        raise ValueError("--sample-dt-fs must be positive")
    if zero_padding < 1:
        raise ValueError("--zero-padding must be at least 1")
    if not 0 <= discard_frames < len(dipole) - 3:
        raise ValueError("Discard interval leaves fewer than four dipole samples")

    signal = np.asarray(dipole[discard_frames:], dtype=float)
    signal -= np.mean(signal, axis=0, keepdims=True)
    acf = autocorrelation_multiple_origins(signal)

    if max_correlation_frames is None:
        correlation_length = len(acf)
    else:
        correlation_length = min(max_correlation_frames, len(acf))
    if correlation_length < 4:
        raise ValueError("Correlation interval contains fewer than four samples")

    acf = acf[:correlation_length]
    times_fs = np.arange(correlation_length) * sample_dt_fs
    window = correlation_window(correlation_length, window_mode)
    windowed_acf = acf * window

    transform_length = zero_padding * correlation_length
    spectral_density = (
        2.0
        * sample_dt_fs
        * np.real(np.fft.rfft(windowed_acf, n=transform_length))
    )
    frequencies_hz = np.fft.rfftfreq(
        transform_length, d=sample_dt_fs * 1.0e-15
    )
    wavenumbers = frequencies_hz / C_CM_PER_S

    # omega^2 differs from wavenumber^2 only by a constant factor.  Since the
    # final spectrum is normalized, wavenumber^2 is numerically convenient.
    ir_raw = wavenumbers**2 * spectral_density
    return {
        "time_fs": times_fs,
        "acf": acf,
        "window": window,
        "windowed_acf": windowed_acf,
        "wavenumber_cm-1": wavenumbers,
        "spectral_density": spectral_density,
        "ir_raw": ir_raw,
    }


def normalize_spectrum(
    wavenumbers: np.ndarray,
    raw_intensity: np.ndarray,
    minimum: float,
    maximum: float,
) -> np.ndarray:
    """Clip finite-sampling negatives and normalize in the requested range."""
    if maximum <= minimum:
        raise ValueError("--max-wavenumber must exceed --min-wavenumber")
    selected = (wavenumbers >= minimum) & (wavenumbers <= maximum)
    if not np.any(selected):
        raise ValueError("Requested wavenumber interval contains no FFT points")

    nonnegative = np.maximum(raw_intensity, 0.0)
    scale = float(np.max(nonnegative[selected]))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(
            "The spectrum has no positive intensity in the requested interval"
        )
    return nonnegative / scale


def find_peaks(
    x: np.ndarray,
    y: np.ndarray,
    minimum: float,
    maximum: float,
    threshold: float,
    min_separation: float,
    maximum_peaks: int,
) -> list[tuple[float, float]]:
    """Locate the strongest separated local maxima without SciPy."""
    local = np.flatnonzero(
        (y[1:-1] > y[:-2])
        & (y[1:-1] >= y[2:])
        & (y[1:-1] >= threshold)
    ) + 1
    local = local[(x[local] >= minimum) & (x[local] <= maximum)]
    ranked = local[np.argsort(y[local])[::-1]]

    selected: list[int] = []
    for index in ranked:
        if all(abs(x[index] - x[other]) >= min_separation for other in selected):
            selected.append(int(index))
        if len(selected) == maximum_peaks:
            break
    return sorted(((float(x[i]), float(y[i])) for i in selected), key=lambda p: p[0])


def save_outputs(
    output_prefix: Path,
    result: dict[str, np.ndarray],
    normalized: np.ndarray,
    minimum: float,
    maximum: float,
) -> tuple[Path, Path, Path | None]:
    """Write the spectrum, correlation function, and plot."""
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    spectrum_path = output_prefix.with_suffix(".csv")
    acf_path = output_prefix.with_name(output_prefix.name + "_acf").with_suffix(
        ".csv"
    )
    plot_path = output_prefix.with_suffix(".png")

    np.savetxt(
        spectrum_path,
        np.column_stack(
            (
                result["wavenumber_cm-1"],
                result["spectral_density"],
                result["ir_raw"],
                normalized,
            )
        ),
        delimiter=",",
        header=(
            "wavenumber_cm-1,spectral_density_arb,"
            "ir_intensity_raw_arb,ir_intensity_normalized"
        ),
        comments="",
    )
    np.savetxt(
        acf_path,
        np.column_stack(
            (
                result["time_fs"],
                result["acf"],
                result["window"],
                result["windowed_acf"],
            )
        ),
        delimiter=",",
        header="time_fs,dipole_acf,window,windowed_dipole_acf",
        comments="",
    )

    written_plot: Path | None = plot_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "WARNING: matplotlib is unavailable; CSV files were written but "
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
            color="#174a7e",
            linewidth=1.5,
        )
        ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
        ax.set_ylabel("Relative IR intensity")
        ax.set_xlim(minimum, maximum)
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=220)
        plt.close(fig)

    return spectrum_path, acf_path, written_plot


def run_self_test() -> int:
    """Verify that a synthetic 3650 cm-1 oscillation is recovered."""
    target = 3650.0
    sample_dt_fs = 0.25
    times_fs = np.arange(0.0, 5000.0, sample_dt_fs)
    frequency_per_fs = target * C_CM_PER_S * 1.0e-15
    dipole = np.sin(2.0 * np.pi * frequency_per_fs * times_fs)[:, None]
    result = calculate_ir(
        dipole=dipole,
        sample_dt_fs=sample_dt_fs,
        discard_frames=0,
        max_correlation_frames=None,
        window_mode="hann",
        zero_padding=4,
    )
    normalized = normalize_spectrum(
        result["wavenumber_cm-1"], result["ir_raw"], 3000.0, 4200.0
    )
    interval = (
        (result["wavenumber_cm-1"] >= 3000.0)
        & (result["wavenumber_cm-1"] <= 4200.0)
    )
    recovered = float(
        result["wavenumber_cm-1"][interval][np.argmax(normalized[interval])]
    )
    if abs(recovered - target) > 3.0:
        raise AssertionError(
            f"Self-test peak is {recovered:.3f} cm^-1, expected {target:.3f}"
        )
    print(
        f"Self-test passed: target={target:.3f} cm^-1, "
        f"recovered={recovered:.3f} cm^-1"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate a relative IR spectrum from one or more corresponding "
            "i-PI dipole extras files."
        )
    )
    parser.add_argument(
        "dipole_files",
        nargs="*",
        type=Path,
        help=(
            "i-PI dipole extras file(s). Multiple files are interpreted as "
            "beads and averaged before calculating the correlation."
        ),
    )
    parser.add_argument(
        "--extra-key",
        default=None,
        help=(
            "JSON extras field to read. Usually polarized_dipole_au or "
            "dipole_au. It is detected automatically when the file contains "
            "only one non-step field."
        ),
    )
    parser.add_argument(
        "--component",
        choices=["auto", "polarized", "x", "y", "z", "isotropic"],
        default="auto",
        help=(
            "Component used for vector dipoles. 'auto' means polarized for "
            "scalar data and isotropic for three-component data."
        ),
    )
    parser.add_argument(
        "--sample-dt-fs",
        type=float,
        default=0.25,
        help="Time between saved dipoles in fs (default: 0.25).",
    )
    parser.add_argument(
        "--discard-fs",
        type=float,
        default=0.0,
        help="Initial trajectory interval to discard in fs (default: 0).",
    )
    parser.add_argument(
        "--max-correlation-fs",
        type=float,
        help=(
            "Maximum correlation time in fs. By default the complete retained "
            "trajectory is used."
        ),
    )
    parser.add_argument(
        "--window",
        choices=["hann", "blackman", "none"],
        default="hann",
        help="Window applied to the autocorrelation (default: hann).",
    )
    parser.add_argument(
        "--zero-padding",
        type=int,
        default=4,
        help="FFT zero-padding factor; this smooths but does not add resolution.",
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
        "--peak-threshold",
        type=float,
        default=0.05,
        help="Minimum normalized intensity reported as a peak.",
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
        default=Path("ir_spectrum"),
        help="Output path prefix (default: ir_spectrum).",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run a synthetic-frequency test without reading an i-PI file.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.self_test:
        return run_self_test()
    if not args.dipole_files:
        raise ValueError("Provide at least one i-PI dipole extras file")
    if args.discard_fs < 0.0:
        raise ValueError("--discard-fs cannot be negative")
    if args.max_correlation_fs is not None and args.max_correlation_fs <= 0.0:
        raise ValueError("--max-correlation-fs must be positive")

    dipole, steps, selected_key = load_bead_average(
        args.dipole_files, args.extra_key, args.component
    )
    discard_frames = int(round(args.discard_fs / args.sample_dt_fs))
    maximum_frames = (
        None
        if args.max_correlation_fs is None
        else int(round(args.max_correlation_fs / args.sample_dt_fs)) + 1
    )
    result = calculate_ir(
        dipole=dipole,
        sample_dt_fs=args.sample_dt_fs,
        discard_frames=discard_frames,
        max_correlation_frames=maximum_frames,
        window_mode=args.window,
        zero_padding=args.zero_padding,
    )
    normalized = normalize_spectrum(
        result["wavenumber_cm-1"],
        result["ir_raw"],
        args.min_wavenumber,
        args.max_wavenumber,
    )
    spectrum_path, acf_path, plot_path = save_outputs(
        args.output_prefix,
        result,
        normalized,
        args.min_wavenumber,
        args.max_wavenumber,
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

    step_stride = "unavailable"
    if steps is not None and len(steps) > 1:
        differences = np.diff(steps)
        if np.allclose(differences, differences[0]):
            step_stride = f"{differences[0]:g}"

    retained_frames = len(dipole) - discard_frames
    duration_fs = max(retained_frames - 1, 0) * args.sample_dt_fs
    correlation_duration_fs = result["time_fs"][-1]
    nominal_resolution = (
        1.0
        / (C_CM_PER_S * correlation_duration_fs * 1.0e-15)
        if correlation_duration_fs > 0.0
        else float("inf")
    )
    print(f"Extras field: {selected_key}")
    print(f"Bead files averaged: {len(args.dipole_files)}")
    print(f"Frames retained: {retained_frames}")
    print(f"Saved-step stride reported by i-PI: {step_stride}")
    print(f"Sample interval: {args.sample_dt_fs:.8g} fs")
    print(f"Retained trajectory duration: {duration_fs:.3f} fs")
    print(f"Correlation duration: {correlation_duration_fs:.3f} fs")
    print(f"Nominal Fourier resolution: {nominal_resolution:.3f} cm^-1")
    print(f"Spectrum: {spectrum_path}")
    print(f"Autocorrelation: {acf_path}")
    if plot_path is not None:
        print(f"Plot: {plot_path}")
    if peaks:
        print("Strongest separated peaks:")
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
