#!/usr/bin/env python3
"""Extract spectral maxima and make a cavity peak-position diagram.

Example
-------
python plot_peak_position_diagram.py \
  --spectrum 0.0186582=Freq0.0186582/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187038=Freq0.0187038/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187505=Freq0.0187505/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187949=Freq0.0187949/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0188404=Freq0.0188404/cavity_IR_lab_z_spectrum.csv \
  --window 4050 4270 --bare-mode 4115 \
  --prominence 0.03 --distance-cm 12 --max-peaks 4 \
  --output water_peak_position.png

Each --spectrum argument is OMEGA_C_AU=CSV_FILE.  The reader accepts comma-,
semicolon-, or whitespace-separated text.  By default it uses the first two
numeric columns as wavenumber and intensity; use --x-col/--y-col if needed.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks


AU_TO_CM = 219474.6313705


def parse_spectrum_spec(text: str) -> tuple[float, Path]:
    try:
        omega_text, filename = text.split("=", 1)
        omega_au = float(omega_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected OMEGA_AU=FILE, received {text!r}"
        ) from exc
    path = Path(filename).expanduser()
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Spectrum file not found: {path}")
    return omega_au, path


def load_spectrum(path: Path, x_col: int | None, y_col: int | None):
    # Header lines become NaN and are discarded after numeric conversion.
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
            f"available numeric columns are {usable}"
        )

    x = numeric[ix].to_numpy(float)
    y = numeric[iy].to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    order = np.argsort(x)
    x, y = x[order], y[order]
    if x.size < 3:
        raise ValueError(f"Too few usable rows in {path}")
    return x, y


def detect_peaks(
    x: np.ndarray,
    y: np.ndarray,
    window: tuple[float, float],
    prominence_fraction: float,
    distance_cm: float,
    max_peaks: int,
):
    select = (x >= window[0]) & (x <= window[1])
    xw, yw = x[select], y[select]
    if xw.size < 3:
        raise ValueError(f"No usable spectrum points inside {window}")

    baseline = np.nanmin(yw)
    scale = np.nanmax(yw) - baseline
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("The selected spectral window has no positive range")
    yn = (yw - baseline) / scale

    dx = float(np.median(np.diff(xw)))
    distance_points = max(1, int(round(distance_cm / abs(dx))))
    indices, properties = find_peaks(
        yn,
        prominence=prominence_fraction,
        distance=distance_points,
    )
    if indices.size == 0:
        return xw, yn, np.array([], dtype=int), np.array([], dtype=float)

    prominences = properties["prominences"]
    if max_peaks > 0 and indices.size > max_peaks:
        keep = np.argsort(prominences)[-max_peaks:]
        indices, prominences = indices[keep], prominences[keep]

    order = np.argsort(xw[indices])
    return xw, yn, indices[order], prominences[order]


def build_parser():
    parser = argparse.ArgumentParser(
        description="Extract maxima from cavity spectra and plot peak position vs cavity frequency."
    )
    parser.add_argument(
        "--spectrum",
        action="append",
        required=True,
        type=parse_spectrum_spec,
        metavar="OMEGA_AU=FILE",
        help="Repeat once for every cavity spectrum.",
    )
    parser.add_argument("--window", nargs=2, type=float, default=(4050.0, 4270.0))
    parser.add_argument(
        "--prominence",
        type=float,
        default=0.03,
        help="Minimum peak prominence as a fraction of that spectrum's window maximum.",
    )
    parser.add_argument("--distance-cm", type=float, default=12.0)
    parser.add_argument(
        "--max-peaks",
        type=int,
        default=4,
        help="Keep this many most-prominent peaks per spectrum; 0 keeps all.",
    )
    parser.add_argument(
        "--bare-mode",
        action="append",
        type=float,
        default=[],
        metavar="CM-1",
        help="Bare molecular peak to draw horizontally; may be repeated.",
    )
    parser.add_argument("--x-col", type=int, default=None, help="Zero-based frequency column.")
    parser.add_argument("--y-col", type=int, default=None, help="Zero-based intensity column.")
    parser.add_argument("--output", default="peak_position_diagram.png")
    parser.add_argument(
        "--peaks-csv",
        default=None,
        help="Output table; default is OUTPUT with suffix _peaks.csv.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def main():
    args = build_parser().parse_args()
    if not 0.0 < args.prominence < 1.0:
        raise SystemExit("--prominence must lie between 0 and 1")
    if args.window[0] >= args.window[1]:
        raise SystemExit("--window requires LOW HIGH")

    spectra = sorted(args.spectrum, key=lambda item: item[0])
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(spectra)))
    results = []
    processed = []

    for (omega_au, path), color in zip(spectra, colors):
        x, y = load_spectrum(path, args.x_col, args.y_col)
        xw, yn, indices, prominences = detect_peaks(
            x,
            y,
            tuple(args.window),
            args.prominence,
            args.distance_cm,
            args.max_peaks,
        )
        cavity_cm = omega_au * AU_TO_CM
        for rank, (idx, prom) in enumerate(zip(indices, prominences), start=1):
            results.append(
                {
                    "omega_c_au": omega_au,
                    "omega_c_cm-1": cavity_cm,
                    "peak_rank_by_frequency": rank,
                    "peak_cm-1": xw[idx],
                    "relative_intensity": yn[idx],
                    "relative_prominence": prom,
                    "source_file": str(path),
                }
            )
        processed.append((omega_au, cavity_cm, xw, yn, indices, color))

    output = Path(args.output)
    peaks_csv = (
        Path(args.peaks_csv)
        if args.peaks_csv
        else output.with_name(output.stem + "_peaks.csv")
    )

    with peaks_csv.open("w", newline="") as handle:
        fieldnames = [
            "omega_c_au",
            "omega_c_cm-1",
            "peak_rank_by_frequency",
            "peak_cm-1",
            "relative_intensity",
            "relative_prominence",
            "source_file",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    fig, (ax_spectra, ax_peaks) = plt.subplots(
        1, 2, figsize=(12.0, 5.2), constrained_layout=True
    )

    # Validation panel: normalized traces plus the exact selected maxima.
    vertical_step = 1.12
    for row, (omega_au, cavity_cm, xw, yn, indices, color) in enumerate(processed):
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
            rf"$\omega_c={omega_au:.7f}$ a.u.",
            color=color,
            va="center",
            fontsize=8,
        )
    ax_spectra.set_xlim(args.window[0], args.window[1] + 62)
    ax_spectra.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    ax_spectra.set_ylabel("Normalized intensity + offset")
    ax_spectra.set_title("Detected maxima (verify this panel)")
    ax_spectra.set_yticks([])

    # Peak-position diagram. Marker size and color both represent IR intensity.
    if results:
        cavity_values = np.array([row["omega_c_cm-1"] for row in results])
        peak_values = np.array([row["peak_cm-1"] for row in results])
        intensities = np.array([row["relative_intensity"] for row in results])
        scatter = ax_peaks.scatter(
            cavity_values,
            peak_values,
            c=intensities,
            s=35 + 115 * intensities,
            cmap="magma",
            vmin=0,
            vmax=1,
            edgecolor="black",
            linewidth=0.45,
            zorder=3,
        )
        colorbar = fig.colorbar(scatter, ax=ax_peaks, pad=0.02)
        colorbar.set_label("Relative IR intensity")

    cavity_scan = np.array([row[1] for row in processed])
    x_padding = max(12.0, 0.15 * np.ptp(cavity_scan))
    xx = np.linspace(cavity_scan.min() - x_padding, cavity_scan.max() + x_padding, 200)
    ax_peaks.plot(xx, xx, ls=":", lw=1.2, color="0.45", label=r"bare photon: $y=x$")
    for i, mode in enumerate(args.bare_mode):
        ax_peaks.axhline(
            mode,
            ls="--",
            lw=1.2,
            color=plt.cm.tab10(i % 10),
            label=rf"bare molecular mode: {mode:.1f} cm$^{{-1}}$",
        )

    ax_peaks.set_xlim(xx.min(), xx.max())
    ax_peaks.set_ylim(args.window)
    ax_peaks.set_xlabel(r"Bare cavity frequency, $\omega_c$ (cm$^{-1}$)")
    ax_peaks.set_ylabel(r"Observed peak position (cm$^{-1}$)")
    ax_peaks.set_title("Peak-position diagram")
    ax_peaks.grid(alpha=0.2)
    ax_peaks.legend(frameon=False, fontsize=8, loc="best")

    fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved figure: {output.resolve()}")
    print(f"Saved peak table: {peaks_csv.resolve()}")
    print(f"Detected {len(results)} peaks across {len(spectra)} spectra.")


if __name__ == "__main__":
    main()
