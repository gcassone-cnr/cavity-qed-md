#!/usr/bin/env python3
"""Make a photon-coordinate q peak-position diagram for a cavity scan.

The script has two mutually exclusive input modes.

1. Read the fictitious photon particle (default label ``L``) directly from
   i-PI extended-XYZ position trajectories. Repeat an omega value to average
   independent replicas before detecting peaks::

    python plot_q_peak_position_diagram.py \
      --trajectory 0.0186582=Freq0.0186582/Sim1/positions_au.xyz \
      --trajectory 0.0186582=Freq0.0186582/Sim2/positions_au.xyz \
      --trajectory 0.0186582=Freq0.0186582/Sim3/positions_au.xyz \
      --trajectory 0.0187038=Freq0.0187038/Sim1/positions_au.xyz \
      --trajectory 0.0187038=Freq0.0187038/Sim2/positions_au.xyz \
      --trajectory 0.0187038=Freq0.0187038/Sim3/positions_au.xyz \
      --dt-fs 0.25 --max-correlation-fs 3000 \
      --broadening-cm 10 --window 4060 4160 \
      --bare-mode 4115.6 --prominence 0.02 \
      --output q_peak_position.png

   ``--dt-fs`` is the time between SAVED XYZ frames, not necessarily the MD
   integration timestep. With ``--q-component auto``, the L-coordinate axis
   having the largest variance is selected independently for every file.

2. Read already calculated q-spectrum files. The first two numeric columns
   are treated as wavenumber and spectral density by default::

    python plot_q_peak_position_diagram.py \
      --spectrum 0.0186582=Freq0.0186582/q_spectrum.csv \
      --spectrum 0.0187038=Freq0.0187038/q_spectrum.csv \
      --spectrum 0.0187505=Freq0.0187505/q_spectrum.csv \
      --spectrum 0.0187949=Freq0.0187949/q_spectrum.csv \
      --spectrum 0.0188404=Freq0.0188404/q_spectrum.csv \
      --window 4060 4160 --bare-mode 4115.6 \
      --prominence 0.02 --output q_peak_position.png

For trajectory input, the plotted spectral density is the cosine transform of
the normalized q autocorrelation function. Peak positions are the same if the
q-velocity spectrum is requested with ``--observable qdot``; only the relative
spectral weights change. The output figure always contains a validation panel,
and the selected peak positions are also written to CSV.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import detrend, find_peaks


AU_TO_CM = 219474.6313705
C_CM_PER_FS = 2.99792458e-5


def parse_input_spec(text: str) -> tuple[float, Path]:
    """Parse OMEGA_C_AU=FILE and verify that FILE exists."""
    try:
        omega_text, filename = text.split("=", 1)
        omega_au = float(omega_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected OMEGA_C_AU=FILE, received {text!r}"
        ) from exc
    path = Path(filename).expanduser()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Input file not found: {path}")
    return omega_au, path


def load_numeric_spectrum(
    path: Path, x_col: int | None, y_col: int | None
) -> tuple[np.ndarray, np.ndarray]:
    """Read a comma-, semicolon-, or whitespace-separated spectrum."""
    raw = pd.read_csv(
        path,
        sep=r"[\s,;]+",
        engine="python",
        comment="#",
        header=None,
    )
    numeric = raw.apply(pd.to_numeric, errors="coerce")
    usable = [int(c) for c in numeric.columns if numeric[c].notna().sum() >= 3]
    if len(usable) < 2:
        raise ValueError(f"Could not find two numeric columns in {path}")

    ix = usable[0] if x_col is None else x_col
    iy = usable[1] if y_col is None else y_col
    if ix not in numeric.columns or iy not in numeric.columns:
        raise ValueError(
            f"Requested columns {ix}, {iy} do not exist in {path}; "
            f"numeric columns are {usable}"
        )

    x = numeric[ix].to_numpy(float)
    y = numeric[iy].to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    order = np.argsort(x)
    x, y = x[order], y[order]
    if x.size < 3 or np.any(np.diff(x) <= 0):
        raise ValueError(f"Spectrum in {path} needs at least 3 distinct x values")
    return x, y


def load_q_from_xyz(
    path: Path, atom_label: str, component: str
) -> tuple[np.ndarray, str]:
    """Extract one coordinate of the fictitious photon atom from extended XYZ."""
    coordinates: list[tuple[float, float, float]] = []
    frame = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        while True:
            first = handle.readline()
            while first and not first.strip():
                first = handle.readline()
            if not first:
                break
            try:
                natoms = int(first.strip().split()[0])
            except (ValueError, IndexError) as exc:
                raise ValueError(
                    f"Invalid XYZ atom-count line in {path}, frame {frame + 1}: "
                    f"{first.rstrip()!r}"
                ) from exc

            comment = handle.readline()
            if not comment:
                raise ValueError(f"Truncated XYZ comment line in {path}, frame {frame + 1}")

            matches: list[tuple[float, float, float]] = []
            for atom_index in range(natoms):
                line = handle.readline()
                if not line:
                    raise ValueError(
                        f"Truncated XYZ atom block in {path}, frame {frame + 1}"
                    )
                fields = line.split()
                if len(fields) < 4:
                    raise ValueError(
                        f"Invalid XYZ atom line in {path}, frame {frame + 1}, "
                        f"atom {atom_index + 1}: {line.rstrip()!r}"
                    )
                if fields[0].casefold() == atom_label.casefold():
                    try:
                        matches.append(tuple(float(value) for value in fields[1:4]))
                    except ValueError as exc:
                        raise ValueError(
                            f"Non-numeric photon coordinates in {path}, frame {frame + 1}"
                        ) from exc

            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly one {atom_label!r} entry in {path}, frame "
                    f"{frame + 1}; found {len(matches)}"
                )
            coordinates.append(matches[0])
            frame += 1

    if len(coordinates) < 8:
        raise ValueError(f"Only {len(coordinates)} XYZ frames found in {path}")

    xyz = np.asarray(coordinates, dtype=float)
    if component == "auto":
        variances = np.var(xyz, axis=0)
        component_index = int(np.argmax(variances))
        selected = "xyz"[component_index]
        if variances[component_index] <= 0:
            raise ValueError(f"Photon coordinates have zero variance in {path}")
    else:
        component_index = "xyz".index(component)
        selected = component

    return xyz[:, component_index], selected


def normalized_autocorrelation(series: np.ndarray, detrend_mode: str) -> np.ndarray:
    """Return the unbiased autocorrelation, normalized so C(0)=1."""
    values = np.asarray(series, dtype=float)
    if detrend_mode == "linear":
        values = detrend(values, type="linear")
    elif detrend_mode == "constant":
        values = values - np.mean(values)

    n = values.size
    nfft = 1 << math.ceil(math.log2(2 * n - 1))
    transform = np.fft.rfft(values, n=nfft)
    acf = np.fft.irfft(transform * np.conjugate(transform), n=nfft)[:n]
    acf /= np.arange(n, 0, -1, dtype=float)
    if not np.isfinite(acf[0]) or acf[0] <= 0:
        raise ValueError("Cannot normalize q autocorrelation because its variance is zero")
    return acf / acf[0]


def correlation_taper(length: int, kind: str) -> np.ndarray:
    """Return a one-sided correlation window equal to one at t=0."""
    if length <= 1 or kind == "none":
        return np.ones(length)
    u = np.linspace(0.0, 1.0, length)
    if kind == "hann":
        return 0.5 * (1.0 + np.cos(np.pi * u))
    if kind == "gaussian":
        # sigma=T/3 gives a weight of exp(-4.5)=0.011 at the cutoff.
        return np.exp(-0.5 * (3.0 * u) ** 2)
    raise ValueError(f"Unknown correlation window: {kind}")


def cosine_spectrum(
    acf: np.ndarray,
    dt_fs: float,
    wavenumbers: np.ndarray,
    observable: str,
    chunk_size: int = 256,
) -> np.ndarray:
    """Cosine-transform a one-sided autocorrelation on a chosen cm^-1 grid."""
    times_fs = np.arange(acf.size, dtype=float) * dt_fs
    integration_weights = np.ones(acf.size)
    if acf.size > 1:
        integration_weights[[0, -1]] = 0.5
    weighted_acf = acf * integration_weights

    spectrum = np.empty(wavenumbers.size, dtype=float)
    phase_scale = 2.0 * np.pi * C_CM_PER_FS
    for start in range(0, wavenumbers.size, chunk_size):
        stop = min(start + chunk_size, wavenumbers.size)
        phases = phase_scale * np.outer(wavenumbers[start:stop], times_fs)
        spectrum[start:stop] = 2.0 * dt_fs * (np.cos(phases) @ weighted_acf)

    if observable == "qdot":
        # S_qdot is proportional to omega^2 S_q. Scaling by max(omega)^2
        # avoids unnecessarily large numbers and does not move the peaks.
        reference = max(float(np.max(wavenumbers)), 1.0)
        spectrum *= (wavenumbers / reference) ** 2
    return spectrum


def trajectory_group_spectrum(
    paths: list[Path], args: argparse.Namespace
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Average normalized replica ACFs and calculate one q spectrum."""
    correlations: list[np.ndarray] = []
    component_messages: list[str] = []
    max_points = int(math.floor(args.max_correlation_fs / args.dt_fs)) + 1

    for path in paths:
        q, selected_component = load_q_from_xyz(
            path, args.atom_label, args.q_component
        )
        acf = normalized_autocorrelation(q, args.detrend)
        usable = min(acf.size, max_points)
        if usable < 8:
            raise ValueError(
                f"Correlation interval leaves only {usable} points for {path}"
            )
        correlations.append(acf[:usable])
        component_messages.append(
            f"{path}: selected {args.atom_label}.{selected_component}; "
            f"{q.size} frames"
        )

    common_length = min(item.size for item in correlations)
    averaged_acf = np.mean(
        np.stack([item[:common_length] for item in correlations], axis=0), axis=0
    )
    averaged_acf *= correlation_taper(common_length, args.correlation_window)

    lo, hi = args.spectrum_range
    wavenumbers = np.arange(lo, hi + 0.5 * args.grid_cm, args.grid_cm)
    spectrum = cosine_spectrum(
        averaged_acf, args.dt_fs, wavenumbers, args.observable
    )

    if args.broadening_cm > 0:
        sigma_cm = args.broadening_cm / 2.354820045
        spectrum = gaussian_filter1d(
            spectrum, sigma=sigma_cm / args.grid_cm, mode="nearest"
        )
    return wavenumbers, spectrum, component_messages


def average_precomputed_spectra(
    paths: list[Path], args: argparse.Namespace
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate and equally average one or more precomputed replica spectra."""
    loaded = [load_numeric_spectrum(path, args.x_col, args.y_col) for path in paths]
    lo = max(float(x[0]) for x, _ in loaded)
    hi = min(float(x[-1]) for x, _ in loaded)
    if lo >= hi:
        raise ValueError("Precomputed replica spectra have no overlapping x range")

    dx = max(float(np.median(np.diff(x))) for x, _ in loaded)
    grid = np.arange(lo, hi + 0.25 * dx, dx)
    replicas = []
    for x, y in loaded:
        interpolated = np.interp(grid, x, y)
        selected = (grid >= args.window[0]) & (grid <= args.window[1])
        if np.count_nonzero(selected) < 3:
            raise ValueError("Requested peak window is absent from a q spectrum")
        baseline = float(np.min(interpolated[selected]))
        scale = float(np.max(interpolated[selected]) - baseline)
        if scale <= 0:
            raise ValueError("A precomputed q spectrum has no range in --window")
        replicas.append((interpolated - baseline) / scale)
    return grid, np.mean(np.stack(replicas, axis=0), axis=0)


def normalized_window(
    x: np.ndarray, y: np.ndarray, window: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    selected = (x >= window[0]) & (x <= window[1])
    xw, yw = x[selected], y[selected]
    if xw.size < 3:
        raise ValueError(f"No usable spectrum points inside {window}")
    baseline = float(np.nanmin(yw))
    scale = float(np.nanmax(yw) - baseline)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Selected spectral window has no positive range")
    return xw, (yw - baseline) / scale


def detect_peaks(
    x: np.ndarray,
    y: np.ndarray,
    window: tuple[float, float],
    prominence_fraction: float,
    distance_cm: float,
    max_peaks: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xw, yn = normalized_window(x, y, window)
    dx = float(np.median(np.diff(xw)))
    distance_points = max(1, int(round(distance_cm / abs(dx))))
    indices, properties = find_peaks(
        yn, prominence=prominence_fraction, distance=distance_points
    )
    if indices.size == 0:
        return xw, yn, np.array([], dtype=int), np.array([], dtype=float)

    prominences = properties["prominences"]
    if max_peaks > 0 and indices.size > max_peaks:
        keep = np.argsort(prominences)[-max_peaks:]
        indices, prominences = indices[keep], prominences[keep]
    order = np.argsort(xw[indices])
    return xw, yn, indices[order], prominences[order]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Calculate/read photon-coordinate q spectra, extract their maxima, "
            "and plot peak position against bare cavity frequency."
        ),
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--trajectory",
        action="append",
        type=parse_input_spec,
        metavar="OMEGA_AU=XYZ",
        help="i-PI position trajectory; repeat an omega to average replicas.",
    )
    inputs.add_argument(
        "--spectrum",
        action="append",
        type=parse_input_spec,
        metavar="OMEGA_AU=CSV",
        help="Precomputed q spectrum; repeat an omega to average replicas.",
    )

    parser.add_argument(
        "--dt-fs",
        type=float,
        default=None,
        help="Time between saved XYZ frames; required with --trajectory.",
    )
    parser.add_argument("--atom-label", default="L")
    parser.add_argument(
        "--q-component", choices=("auto", "x", "y", "z"), default="auto"
    )
    parser.add_argument(
        "--detrend", choices=("linear", "constant", "none"), default="linear"
    )
    parser.add_argument(
        "--observable",
        choices=("q", "qdot"),
        default="q",
        help="q-coordinate or q-velocity-weighted spectral density.",
    )
    parser.add_argument("--max-correlation-fs", type=float, default=3000.0)
    parser.add_argument(
        "--correlation-window",
        choices=("gaussian", "hann", "none"),
        default="gaussian",
    )
    parser.add_argument(
        "--broadening-cm",
        type=float,
        default=10.0,
        help="Additional Gaussian FWHM applied to trajectory-derived spectra.",
    )
    parser.add_argument(
        "--spectrum-range",
        nargs=2,
        type=float,
        default=(3800.0, 4400.0),
        metavar=("LOW", "HIGH"),
        help="Calculated range for trajectory-derived q spectra.",
    )
    parser.add_argument("--grid-cm", type=float, default=1.0)

    parser.add_argument(
        "--window", nargs=2, type=float, default=(4050.0, 4270.0)
    )
    parser.add_argument("--prominence", type=float, default=0.02)
    parser.add_argument("--distance-cm", type=float, default=8.0)
    parser.add_argument(
        "--max-peaks",
        type=int,
        default=4,
        help="Most-prominent peaks retained per omega; 0 retains all.",
    )
    parser.add_argument(
        "--bare-mode",
        action="append",
        type=float,
        default=[],
        metavar="CM-1",
    )
    parser.add_argument("--x-col", type=int, default=None)
    parser.add_argument("--y-col", type=int, default=None)
    parser.add_argument("--output", default="q_peak_position.png")
    parser.add_argument("--peaks-csv", default=None)
    parser.add_argument(
        "--spectra-dir",
        default=None,
        help="Directory for trajectory-derived q-spectrum CSV files.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.trajectory and (args.dt_fs is None or args.dt_fs <= 0):
        raise SystemExit("--dt-fs must be positive and is required with --trajectory")
    if not 0.0 < args.prominence < 1.0:
        raise SystemExit("--prominence must lie between 0 and 1")
    if args.window[0] >= args.window[1]:
        raise SystemExit("--window requires LOW HIGH")
    if args.spectrum_range[0] >= args.spectrum_range[1]:
        raise SystemExit("--spectrum-range requires LOW HIGH")
    if args.trajectory and (
        args.window[0] < args.spectrum_range[0]
        or args.window[1] > args.spectrum_range[1]
    ):
        raise SystemExit("--window must lie inside --spectrum-range")
    if args.max_correlation_fs <= 0:
        raise SystemExit("--max-correlation-fs must be positive")
    if args.grid_cm <= 0 or args.broadening_cm < 0 or args.distance_cm <= 0:
        raise SystemExit("Grid and distance must be positive; broadening cannot be negative")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)

    source_specs = args.trajectory if args.trajectory else args.spectrum
    grouped: dict[float, list[Path]] = defaultdict(list)
    for omega_au, path in source_specs:
        grouped[omega_au].append(path)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    peaks_csv = (
        Path(args.peaks_csv)
        if args.peaks_csv
        else output.with_name(output.stem + "_peaks.csv")
    )
    peaks_csv.parent.mkdir(parents=True, exist_ok=True)

    spectra_dir: Path | None = None
    if args.trajectory:
        spectra_dir = (
            Path(args.spectra_dir)
            if args.spectra_dir
            else output.with_name(output.stem + "_spectra")
        )
        spectra_dir.mkdir(parents=True, exist_ok=True)

    processed = []
    results: list[dict[str, object]] = []

    for omega_au in sorted(grouped):
        paths = grouped[omega_au]
        if args.trajectory:
            x, y, messages = trajectory_group_spectrum(paths, args)
            for message in messages:
                print(message)
            spectrum_file = spectra_dir / f"q_spectrum_omega_{omega_au:.7f}.csv"
            pd.DataFrame(
                {"wavenumber_cm-1": x, "q_spectral_density": y}
            ).to_csv(spectrum_file, index=False)
        else:
            x, y = average_precomputed_spectra(paths, args)
            spectrum_file = None

        xw, yn, indices, prominences = detect_peaks(
            x,
            y,
            tuple(args.window),
            args.prominence,
            args.distance_cm,
            args.max_peaks,
        )
        cavity_cm = omega_au * AU_TO_CM
        for rank, (idx, prominence) in enumerate(
            zip(indices, prominences), start=1
        ):
            results.append(
                {
                    "omega_c_au": omega_au,
                    "omega_c_cm-1": cavity_cm,
                    "peak_rank_by_frequency": rank,
                    "peak_cm-1": float(xw[idx]),
                    "relative_q_spectral_weight": float(yn[idx]),
                    "relative_prominence": float(prominence),
                    "replicas": len(paths),
                    "source_files": " | ".join(str(path) for path in paths),
                }
            )
        processed.append((omega_au, cavity_cm, xw, yn, indices, len(paths)))

    with peaks_csv.open("w", newline="") as handle:
        fieldnames = [
            "omega_c_au",
            "omega_c_cm-1",
            "peak_rank_by_frequency",
            "peak_cm-1",
            "relative_q_spectral_weight",
            "relative_prominence",
            "replicas",
            "source_files",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(processed)))
    fig, (ax_spectra, ax_peaks) = plt.subplots(
        1, 2, figsize=(12.0, 5.2), constrained_layout=True
    )

    vertical_step = 1.12
    for row, (record, color) in enumerate(zip(processed, colors)):
        omega_au, _cavity_cm, xw, yn, indices, replicas = record
        offset = row * vertical_step
        ax_spectra.plot(xw, yn + offset, color=color, lw=1.5)
        ax_spectra.scatter(
            xw[indices],
            yn[indices] + offset,
            s=35,
            facecolor="white",
            edgecolor=color,
            linewidth=1.5,
            zorder=3,
        )
        ax_spectra.text(
            args.window[1] + 3,
            offset + 0.50,
            rf"$\omega_c={omega_au:.7f}$ a.u. ({replicas} rep.)",
            color=color,
            va="center",
            fontsize=8,
        )
    ax_spectra.set_xlim(args.window[0], args.window[1] + 72)
    ax_spectra.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax_spectra.set_ylabel(r"Normalized $q$ spectral density + offset")
    ax_spectra.set_title(r"Detected maxima in $q$ spectra (verify this panel)")
    ax_spectra.set_yticks([])

    if results:
        cavity_values = np.array([row["omega_c_cm-1"] for row in results])
        peak_values = np.array([row["peak_cm-1"] for row in results])
        weights = np.array(
            [row["relative_q_spectral_weight"] for row in results]
        )
        scatter = ax_peaks.scatter(
            cavity_values,
            peak_values,
            c=weights,
            s=35 + 115 * weights,
            cmap="viridis",
            vmin=0,
            vmax=1,
            edgecolor="black",
            linewidth=0.45,
            zorder=3,
        )
        colorbar = fig.colorbar(scatter, ax=ax_peaks, pad=0.02)
        colorbar.set_label(r"Relative $q$ spectral weight")

    cavity_scan = np.array([row[1] for row in processed])
    x_padding = max(12.0, 0.15 * float(np.ptp(cavity_scan)))
    xx = np.linspace(
        float(cavity_scan.min()) - x_padding,
        float(cavity_scan.max()) + x_padding,
        200,
    )
    ax_peaks.plot(
        xx, xx, ls=":", lw=1.2, color="0.45", label=r"bare photon: $y=x$"
    )
    for index, mode in enumerate(args.bare_mode):
        ax_peaks.axhline(
            mode,
            ls="--",
            lw=1.2,
            color=plt.cm.tab10(index % 10),
            label=rf"bare molecular mode: {mode:.1f} cm$^{{-1}}$",
        )

    ax_peaks.set_xlim(xx.min(), xx.max())
    ax_peaks.set_ylim(args.window)
    ax_peaks.set_xlabel(r"Bare cavity frequency, $\omega_c$ (cm$^{-1}$)")
    ax_peaks.set_ylabel(r"Observed $q$-spectrum peak (cm$^{-1}$)")
    ax_peaks.set_title(r"Photon-coordinate $q$ peak-position diagram")
    ax_peaks.grid(alpha=0.2)
    ax_peaks.legend(frameon=False, fontsize=8, loc="best")

    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved figure: {output.resolve()}")
    print(f"Saved peak table: {peaks_csv.resolve()}")
    if spectra_dir is not None:
        print(f"Saved q spectra: {spectra_dir.resolve()}")
    print(f"Detected {len(results)} peaks across {len(processed)} cavity frequencies.")


if __name__ == "__main__":
    main()
