#!/usr/bin/env python3
"""Average independent IR spectra on a common frequency grid.

The program is intended for independent MD replicas.  It interpolates the
*raw* intensity from every input spectrum onto one common grid, calculates the
mean, sample standard deviation (SD), and standard error of the mean (SEM), and
only then normalizes the averaged result.  An optional Gaussian broadening is
applied to each interpolated replica before calculating the smoothed mean and
uncertainty.

Typical use
-----------
python3 average_ir_spectra.py rep1.csv rep2.csv rep3.csv \
    --range 3500 4300 --normalization-range 3500 4300 \
    --smooth-fwhm-cm 6 --output-prefix ir_cavity_average

By default, a named column containing "raw" and "intensity" is selected.  For
the CSV files produced by ir_spectrum.py this is ir_intensity_raw_arb (column
3).  A 1-based column number can be selected explicitly with, for example,
--intensity-column 3.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


@dataclass
class Spectrum:
    path: Path
    wavenumber: np.ndarray
    intensity: np.ndarray
    wavenumber_label: str
    intensity_label: str


def _split_line(line: str, delimiter: str | None) -> list[str]:
    line = line.strip().lstrip("#").strip()
    if delimiter is None:
        return line.split()
    return next(csv.reader([line], delimiter=delimiter))


def _detect_layout(path: Path) -> tuple[str | None, list[str] | None, int]:
    """Return delimiter, optional header fields, and number of header rows."""
    first = ""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                first = stripped
                break
    if not first:
        raise ValueError(f"{path}: file is empty")

    if "," in first:
        delimiter: str | None = ","
    elif "\t" in first:
        delimiter = "\t"
    else:
        delimiter = None

    tokens = _split_line(first, delimiter)
    try:
        for token in tokens:
            float(token)
        has_header = False
    except ValueError:
        has_header = True

    return delimiter, tokens if has_header else None, 1 if has_header else 0


def _canonical(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def _resolve_column(
    requested: str,
    header: list[str] | None,
    ncols: int,
    kind: str,
) -> int:
    """Resolve a 1-based number, exact/partial label, or an automatic column."""
    if requested.lower() != "auto":
        try:
            index = int(requested) - 1
        except ValueError:
            if header is None:
                raise ValueError(
                    f"Cannot select column {requested!r} by name: file has no header"
                )
            target = _canonical(requested)
            canon = [_canonical(item) for item in header]
            exact = [i for i, name in enumerate(canon) if name == target]
            partial = [i for i, name in enumerate(canon) if target in name]
            matches = exact or partial
            if len(matches) != 1:
                raise ValueError(
                    f"Column name {requested!r} matched {len(matches)} columns: {header}"
                )
            return matches[0]
        if not 0 <= index < ncols:
            raise ValueError(
                f"Requested {kind} column {index + 1}, but file has {ncols} columns"
            )
        return index

    if header is not None:
        canon = [_canonical(item) for item in header]
        if kind == "wavenumber":
            priorities = (
                lambda s: "wavenumber" in s,
                lambda s: s in {"wn", "frequencycm1", "frequency"},
                lambda s: "frequency" in s,
            )
        else:
            priorities = (
                lambda s: "irintensityraw" in s,
                lambda s: "meanraw" in s,
                lambda s: "intensity" in s and "raw" in s,
                lambda s: "intensity" in s and "norm" not in s,
                lambda s: s in {"raw", "signal", "spectrum"},
            )
        for predicate in priorities:
            matches = [i for i, name in enumerate(canon) if predicate(name)]
            if matches:
                return matches[0]

    if kind == "wavenumber":
        return 0
    if ncols >= 3:
        return 2
    if ncols >= 2:
        return 1
    raise ValueError("At least two numeric columns are required")


def _deduplicate_sorted(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    unique_x, inverse, counts = np.unique(x, return_inverse=True, return_counts=True)
    if unique_x.size != x.size:
        sums = np.zeros(unique_x.size, dtype=float)
        np.add.at(sums, inverse, y)
        y = sums / counts
        x = unique_x
    return x, y


def read_spectrum(
    path: Path,
    wavenumber_column: str,
    intensity_column: str,
) -> Spectrum:
    delimiter, header, skip_header = _detect_layout(path)
    data = np.genfromtxt(
        path,
        delimiter=delimiter,
        comments="#",
        skip_header=skip_header,
        dtype=float,
        invalid_raise=False,
    )
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"{path}: could not read at least two numeric columns")

    wn_col = _resolve_column(wavenumber_column, header, data.shape[1], "wavenumber")
    in_col = _resolve_column(intensity_column, header, data.shape[1], "intensity")
    if wn_col == in_col:
        raise ValueError(f"{path}: wavenumber and intensity resolve to the same column")

    finite = np.isfinite(data[:, wn_col]) & np.isfinite(data[:, in_col])
    x, y = _deduplicate_sorted(data[finite, wn_col], data[finite, in_col])
    if x.size < 3:
        raise ValueError(f"{path}: fewer than three valid spectral points")
    if np.any(np.diff(x) <= 0):
        raise ValueError(f"{path}: wavenumbers are not strictly increasing")

    wn_label = header[wn_col].strip() if header else f"column_{wn_col + 1}"
    in_label = header[in_col].strip() if header else f"column_{in_col + 1}"
    return Spectrum(path, x, y, wn_label, in_label)


def gaussian_smooth(y: np.ndarray, spacing: float, fwhm: float) -> np.ndarray:
    if fwhm <= 0:
        return y.copy()
    sigma_points = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)) * spacing)
    if sigma_points < 0.15:
        return y.copy()
    half_width = max(1, int(math.ceil(4.0 * sigma_points)))
    grid = np.arange(-half_width, half_width + 1, dtype=float)
    kernel = np.exp(-0.5 * (grid / sigma_points) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(y, half_width, mode="reflect")
    return np.convolve(padded, kernel, mode="same")[half_width:-half_width]


def _window_mask(x: np.ndarray, bounds: Sequence[float] | None) -> np.ndarray:
    if bounds is None:
        return np.ones(x.size, dtype=bool)
    low, high = bounds
    return (x >= low) & (x <= high)


def _normalization_divisor(
    x: np.ndarray,
    mean_raw: np.ndarray,
    method: str,
    bounds: Sequence[float] | None,
    supplied_value: float | None,
) -> float:
    if supplied_value is not None:
        if not np.isfinite(supplied_value) or supplied_value <= 0:
            raise ValueError("--normalization-value must be positive and finite")
        return supplied_value
    if method == "none":
        return 1.0

    mask = _window_mask(x, bounds)
    if np.count_nonzero(mask) < 2:
        raise ValueError("Normalization range contains fewer than two grid points")
    xw = x[mask]
    yw = mean_raw[mask]
    if method == "max":
        divisor = float(np.max(yw))
    elif method == "area":
        divisor = float(np.trapezoid(np.clip(yw, 0.0, None), xw))
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    if not np.isfinite(divisor) or divisor <= 0:
        raise ValueError(f"Normalization divisor is not positive: {divisor}")
    return divisor


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Average raw IR spectra from independent MD replicas, calculate "
            "SD/SEM, and normalize only after averaging."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("spectra", nargs="+", type=Path, help="Input CSV/DAT files")
    parser.add_argument(
        "-o", "--output-prefix", type=Path, default=Path("ir_average"),
        help="Output path without extension",
    )
    parser.add_argument(
        "--wavenumber-column", default="auto",
        help="1-based column number, header name, or 'auto'",
    )
    parser.add_argument(
        "--intensity-column", default="auto",
        help="1-based raw-intensity column number, header name, or 'auto'",
    )
    parser.add_argument(
        "--range", nargs=2, type=float, metavar=("MIN", "MAX"),
        help="Restrict the output/common grid to this range in cm^-1",
    )
    parser.add_argument(
        "--normalization", choices=("max", "area", "none"), default="max",
        help="Normalization applied after replica averaging",
    )
    parser.add_argument(
        "--normalization-range", nargs=2, type=float, metavar=("MIN", "MAX"),
        help="Range used to determine the normalization divisor",
    )
    parser.add_argument(
        "--normalization-value", type=float,
        help=(
            "Use this common divisor instead of deriving one.  This is useful "
            "when cavity and no-cavity averages must share the no-cavity scale."
        ),
    )
    parser.add_argument(
        "--smooth-fwhm-cm", type=float, default=0.0,
        help="Optional Gaussian FWHM in cm^-1, applied after interpolation",
    )
    parser.add_argument(
        "--no-plot", action="store_true", help="Do not create a PNG plot",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if len(args.spectra) < 2:
        raise ValueError("Provide at least two independent replica spectra")
    if args.range and args.range[0] >= args.range[1]:
        raise ValueError("--range requires MIN < MAX")
    if args.normalization_range and args.normalization_range[0] >= args.normalization_range[1]:
        raise ValueError("--normalization-range requires MIN < MAX")
    if args.smooth_fwhm_cm < 0:
        raise ValueError("--smooth-fwhm-cm cannot be negative")

    spectra = [
        read_spectrum(path, args.wavenumber_column, args.intensity_column)
        for path in args.spectra
    ]

    lower = max(spec.wavenumber[0] for spec in spectra)
    upper = min(spec.wavenumber[-1] for spec in spectra)
    if args.range:
        lower = max(lower, args.range[0])
        upper = min(upper, args.range[1])
    if lower >= upper:
        raise ValueError("The input spectra have no overlapping requested range")

    first = spectra[0]
    grid_mask = (first.wavenumber >= lower) & (first.wavenumber <= upper)
    common_x = first.wavenumber[grid_mask]
    if common_x.size < 3:
        raise ValueError("Common grid has fewer than three points")
    spacing = float(np.median(np.diff(common_x)))

    replica_raw = np.vstack(
        [np.interp(common_x, spec.wavenumber, spec.intensity) for spec in spectra]
    )
    nrep = replica_raw.shape[0]
    mean_raw = np.mean(replica_raw, axis=0)
    sd_raw = np.std(replica_raw, axis=0, ddof=1)
    sem_raw = sd_raw / math.sqrt(nrep)

    replica_smooth = np.vstack(
        [gaussian_smooth(row, spacing, args.smooth_fwhm_cm) for row in replica_raw]
    )
    mean_smooth = np.mean(replica_smooth, axis=0)
    sd_smooth = np.std(replica_smooth, axis=0, ddof=1)
    sem_smooth = sd_smooth / math.sqrt(nrep)

    divisor = _normalization_divisor(
        common_x,
        mean_raw,
        args.normalization,
        args.normalization_range,
        args.normalization_value,
    )

    output = np.column_stack(
        (
            common_x,
            mean_raw,
            sd_raw,
            sem_raw,
            mean_raw / divisor,
            sem_raw / divisor,
            mean_smooth,
            sd_smooth,
            sem_smooth,
            mean_smooth / divisor,
            sem_smooth / divisor,
        )
    )
    columns = [
        "wavenumber_cm-1",
        "mean_raw",
        "sd_raw",
        "sem_raw",
        "mean_normalized",
        "sem_normalized",
        "mean_smoothed_raw",
        "sd_smoothed_raw",
        "sem_smoothed_raw",
        "mean_smoothed_normalized",
        "sem_smoothed_normalized",
    ]

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    dat_path = prefix.with_suffix(".dat")
    csv_path = prefix.with_suffix(".csv")
    report_path = prefix.with_suffix(".txt")
    np.savetxt(dat_path, output, fmt="%.10e", header=" ".join(columns))
    np.savetxt(
        csv_path,
        output,
        fmt="%.10e",
        delimiter=",",
        header=",".join(columns),
        comments="",
    )

    peak_source = mean_smooth if args.smooth_fwhm_cm > 0 else mean_raw
    norm_mask = _window_mask(common_x, args.normalization_range)
    peak_local = int(np.argmax(peak_source[norm_mask]))
    peak_index = np.flatnonzero(norm_mask)[peak_local]
    report_lines = [
        "IR replica-average report",
        f"replicas: {nrep}",
        f"common range (cm^-1): {common_x[0]:.6f} to {common_x[-1]:.6f}",
        f"median grid spacing (cm^-1): {spacing:.8g}",
        f"normalization method: {args.normalization}",
        f"normalization divisor: {divisor:.12e}",
        f"Gaussian FWHM (cm^-1): {args.smooth_fwhm_cm:.6g}",
        f"strongest plotted peak (cm^-1): {common_x[peak_index]:.6f}",
        "",
        "Input columns:",
    ]
    for spec in spectra:
        report_lines.append(
            f"  {spec.path}: {spec.wavenumber_label!r}, {spec.intensity_label!r}"
        )
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    png_path = prefix.with_suffix(".png")
    plot_written = False
    if not args.no_plot:
        try:
            os.environ.setdefault(
                "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "average_ir_mpl_cache")
            )
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            y = mean_smooth / divisor if args.smooth_fwhm_cm > 0 else mean_raw / divisor
            yerr = sem_smooth / divisor if args.smooth_fwhm_cm > 0 else sem_raw / divisor
            fig, ax = plt.subplots(figsize=(9.0, 5.4))
            ax.fill_between(
                common_x, y - yerr, y + yerr,
                color="#4477AA", alpha=0.25, linewidth=0, label="SEM",
            )
            ax.plot(common_x, y, color="#114477", linewidth=1.8, label=f"Mean (n={nrep})")
            ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
            ylabel = "Normalized IR intensity" if args.normalization != "none" else "Mean raw IR intensity"
            ax.set_ylabel(ylabel)
            ax.set_xlim(common_x[0], common_x[-1])
            ax.grid(alpha=0.22)
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(png_path, dpi=220)
            plt.close(fig)
            plot_written = True
        except ImportError:
            print("matplotlib is unavailable; skipping PNG output", file=sys.stderr)

    print(f"Replicas averaged: {nrep}")
    print(f"Normalization divisor: {divisor:.12e}")
    print(f"Written: {dat_path}")
    print(f"Written: {csv_path}")
    print(f"Written: {report_path}")
    if plot_written:
        print(f"Written: {png_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
