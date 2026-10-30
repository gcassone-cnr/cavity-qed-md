#!/usr/bin/env python3
"""Compare cavity and no-cavity IR spectra and assign a candidate Rabi gap.

The default normalization divides BOTH spectra by one number: the maximum of
the smoothed no-cavity spectrum in the normalization window.  This preserves
their relative intensity scale.  An independently normalized panel is also
drawn, but it is intended only for comparing peak positions and line shapes.

Example
-------
python3 compare_ir_spectra.py ir_cavity.csv ir_no_cavity.csv \
    --omega-au 0.0178191 --trajectory-ps 5 --output-prefix ir_compare
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
except ImportError as exc:  # pragma: no cover - dependency check
    raise SystemExit(
        "This script requires scipy. Install numpy, scipy, and matplotlib in "
        "the Python environment used to run it."
    ) from exc


AU_TO_CM = 219474.6313632
PS_TO_RESOLUTION_CM = 33.3564095198


@dataclass
class Peak:
    wavenumber_cm: float
    height: float
    prominence: float


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Smooth and compare cavity/no-cavity IR spectra using a shared "
            "normalization, then identify a provisional LP/UP pair."
        )
    )
    p.add_argument("cavity", type=Path, help="CSV spectrum with lambda > 0")
    p.add_argument("outside", type=Path, help="CSV spectrum with lambda = 0")
    freq = p.add_mutually_exclusive_group(required=True)
    freq.add_argument("--omega-au", type=float, help="Cavity frequency in atomic units")
    freq.add_argument("--omega-cm", type=float, help="Cavity frequency in cm^-1")
    p.add_argument(
        "--intensity-column",
        default="auto",
        help=(
            "Intensity header to use. 'auto' prefers the raw IR-intensity "
            "column and deliberately avoids an independently normalized column."
        ),
    )
    p.add_argument(
        "--analysis-range",
        nargs=2,
        type=float,
        metavar=("MIN_CM", "MAX_CM"),
        default=(3500.0, 4300.0),
    )
    p.add_argument(
        "--normalization-range",
        nargs=2,
        type=float,
        metavar=("MIN_CM", "MAX_CM"),
        default=None,
        help="Default: same as --analysis-range",
    )
    p.add_argument(
        "--normalization",
        choices=("reference-max", "reference-area", "independent-max", "none"),
        default="reference-max",
        help=(
            "reference-max (recommended) divides both spectra by the no-cavity "
            "maximum; reference-area uses its integrated area; independent-max "
            "is only suitable for line-shape/peak-position comparison."
        ),
    )
    p.add_argument(
        "--smooth-fwhm-cm",
        type=float,
        default=8.0,
        help="Gaussian smoothing FWHM in cm^-1 (default: 8)",
    )
    p.add_argument(
        "--trajectory-ps",
        type=float,
        default=None,
        help="Trajectory length; used only to report the physical FFT resolution",
    )
    p.add_argument(
        "--prominence-fraction",
        type=float,
        default=0.02,
        help="Peak prominence threshold relative to the spectrum maximum",
    )
    p.add_argument(
        "--min-peak-distance-cm",
        type=float,
        default=8.0,
        help="Minimum separation between detected local maxima",
    )
    p.add_argument(
        "--bare-search-halfwidth-cm",
        type=float,
        default=120.0,
        help="Find the bare molecular peak within this distance of omega_c",
    )
    p.add_argument(
        "--rabi-window-cm",
        type=float,
        default=200.0,
        help="Search for LP and UP within +/- this distance of their midpoint",
    )
    p.add_argument(
        "--min-rabi-gap-cm",
        type=float,
        default=30.0,
        help="Reject LP/UP pairs narrower than this value (default: 30 cm^-1)",
    )
    p.add_argument(
        "--max-rabi-gap-cm",
        type=float,
        default=200.0,
        help="Reject LP/UP pairs wider than this value (default: 200 cm^-1)",
    )
    p.add_argument(
        "--center-exclusion-cm",
        type=float,
        default=None,
        help=(
            "Minimum distance of LP/UP from the uncoupled frequencies. "
            "Default: max(smoothing FWHM, reported FFT resolution)."
        ),
    )
    p.add_argument(
        "--midpoint-tolerance-cm",
        type=float,
        default=8.0,
        help="Scale used to reward LP/UP pairs centered on (omega_c + omega_v)/2",
    )
    p.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("ir_comparison"),
        help="Prefix for .dat, .csv, .json, .txt, and .png outputs",
    )
    p.add_argument("--no-plot", action="store_true", help="Do not create PNG output")
    return p.parse_args()


def clean_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def choose_column(names: tuple[str, ...], requested: str) -> str:
    if requested != "auto":
        if requested in names:
            return requested
        requested_clean = clean_name(requested)
        for name in names:
            if clean_name(name) == requested_clean:
                return name
        raise ValueError(
            f"Intensity column '{requested}' not found. Available columns: {', '.join(names)}"
        )

    ranked: list[tuple[int, str]] = []
    for name in names:
        n = clean_name(name)
        score = 0
        if "intensity" in n:
            score += 5
        if "ir" in n:
            score += 2
        if "raw" in n:
            score += 10
        if "normalized" in n or "normalised" in n or n.endswith("norm"):
            score -= 20
        if "wavenumber" in n or "frequency" in n:
            score -= 30
        ranked.append((score, name))
    score, name = max(ranked)
    if score <= 0:
        raise ValueError(
            "Could not identify an intensity column automatically. Use --intensity-column."
        )
    return name


def choose_wavenumber_column(names: tuple[str, ...]) -> str:
    for name in names:
        n = clean_name(name)
        if "wavenumber" in n or n in {"cm1", "frequencycm1", "wncm1"}:
            return name
    return names[0]


def load_spectrum(path: Path, intensity_column: str) -> tuple[np.ndarray, np.ndarray, str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    first_line = path.open("r", encoding="utf-8", errors="replace").readline()
    delimiter = "," if "," in first_line else None
    data = np.genfromtxt(path, delimiter=delimiter, names=True, encoding="utf-8")
    if data.dtype.names is None or len(data.dtype.names) < 2:
        raise ValueError(f"{path}: expected a header and at least two numeric columns")
    names = data.dtype.names
    x_name = choose_wavenumber_column(names)
    y_name = choose_column(names, intensity_column)
    x = np.asarray(data[x_name], dtype=float)
    y = np.asarray(data[y_name], dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    order = np.argsort(x)
    x, y = x[order], y[order]
    unique = np.concatenate(([True], np.diff(x) > 0))
    return x[unique], y[unique], y_name


def make_common_grid(
    xc: np.ndarray,
    yc: np.ndarray,
    xo: np.ndarray,
    yo: np.ndarray,
    bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lo = max(bounds[0], float(xc[0]), float(xo[0]))
    hi = min(bounds[1], float(xc[-1]), float(xo[-1]))
    if hi <= lo:
        raise ValueError("The spectra do not overlap in the requested analysis range")
    dxc = float(np.median(np.diff(xc)))
    dxo = float(np.median(np.diff(xo)))
    dx = max(dxc, dxo)  # do not pretend to have a finer grid than either input
    x = np.arange(lo, hi + 0.25 * dx, dx)
    return x, np.interp(x, xc, yc), np.interp(x, xo, yo)


def smooth(y: np.ndarray, dx: float, fwhm_cm: float) -> np.ndarray:
    if fwhm_cm <= 0:
        return y.copy()
    sigma_points = (fwhm_cm / 2.354820045) / dx
    return gaussian_filter1d(y, sigma=max(sigma_points, 1e-12), mode="nearest")


def window_mask(x: np.ndarray, bounds: tuple[float, float]) -> np.ndarray:
    mask = (x >= bounds[0]) & (x <= bounds[1])
    if np.count_nonzero(mask) < 3:
        raise ValueError(f"Too few points in window {bounds[0]}--{bounds[1]} cm^-1")
    return mask


def normalize(
    x: np.ndarray,
    cavity: np.ndarray,
    outside: np.ndarray,
    bounds: tuple[float, float],
    method: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    mask = window_mask(x, bounds)
    c_pos = np.maximum(cavity, 0.0)
    o_pos = np.maximum(outside, 0.0)
    if method == "reference-max":
        scale = float(np.max(o_pos[mask]))
        scales = {"cavity": scale, "outside": scale}
    elif method == "reference-area":
        scale = float(np.trapezoid(o_pos[mask], x[mask]))
        scales = {"cavity": scale, "outside": scale}
    elif method == "independent-max":
        scales = {
            "cavity": float(np.max(c_pos[mask])),
            "outside": float(np.max(o_pos[mask])),
        }
    else:
        scales = {"cavity": 1.0, "outside": 1.0}
    if scales["cavity"] <= 0 or scales["outside"] <= 0:
        raise ValueError("Cannot normalize: selected reference intensity is not positive")
    return cavity / scales["cavity"], outside / scales["outside"], scales


def detect_peaks(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    prominence_fraction: float,
    min_distance_cm: float,
) -> list[Peak]:
    idx = np.flatnonzero(mask)
    local = np.maximum(y[mask], 0.0)
    maximum = float(np.max(local))
    if maximum <= 0:
        return []
    dx = float(np.median(np.diff(x)))
    distance = max(1, int(round(min_distance_cm / dx)))
    found, props = find_peaks(
        local,
        prominence=prominence_fraction * maximum,
        distance=distance,
    )
    peaks = [
        Peak(
            wavenumber_cm=float(x[idx[i]]),
            height=float(y[idx[i]]),
            prominence=float(props["prominences"][j]),
        )
        for j, i in enumerate(found)
    ]
    return sorted(peaks, key=lambda p: p.wavenumber_cm)


def nearest_peak(peaks: list[Peak], target: float) -> tuple[Peak | None, float | None]:
    if not peaks:
        return None, None
    peak = min(peaks, key=lambda p: abs(p.wavenumber_cm - target))
    return peak, abs(peak.wavenumber_cm - target)


def choose_bare_peak(
    peaks: list[Peak], omega_cm: float, halfwidth: float
) -> Peak | None:
    candidates = [p for p in peaks if abs(p.wavenumber_cm - omega_cm) <= halfwidth]
    if not candidates:
        return None
    # The user-supplied cavity frequency anchors the transition assignment.
    # Prominence is only a deterministic tie-breaker on the discrete FFT grid;
    # otherwise a strong neighbouring rovibrational line can hijack the fit.
    return min(candidates, key=lambda p: (abs(p.wavenumber_cm - omega_cm), -p.prominence))


def choose_polariton_pair(
    peaks: list[Peak],
    omega_cm: float,
    bare_cm: float,
    exclusion_cm: float,
    rabi_window_cm: float,
    midpoint_tolerance_cm: float,
    min_gap_cm: float,
    max_gap_cm: float,
) -> tuple[Peak, Peak, float, list[dict[str, float]]] | None:
    expected_mid = 0.5 * (omega_cm + bare_cm)
    lower_limit = min(omega_cm, bare_cm) - exclusion_cm
    upper_limit = max(omega_cm, bare_cm) + exclusion_cm
    lower = [
        p
        for p in peaks
        if expected_mid - rabi_window_cm <= p.wavenumber_cm < lower_limit
    ]
    upper = [
        p
        for p in peaks
        if upper_limit < p.wavenumber_cm <= expected_mid + rabi_window_cm
    ]
    if not lower or not upper:
        return None
    max_prom = max(p.prominence for p in peaks)
    ranked: list[tuple[float, Peak, Peak, float]] = []
    for lp in lower:
        for up in upper:
            gap = up.wavenumber_cm - lp.wavenumber_cm
            if gap < min_gap_cm or gap > max_gap_cm:
                continue
            midpoint = 0.5 * (lp.wavenumber_cm + up.wavenumber_cm)
            midpoint_error = midpoint - expected_mid
            prominence_term = math.sqrt(lp.prominence * up.prominence) / max_prom
            center_term = math.exp(
                -0.5 * (midpoint_error / max(midpoint_tolerance_cm, 1e-12)) ** 2
            )
            # Centering is the defining constraint; prominence prevents selecting
            # a symmetric pair made only of tiny numerical ripples.
            score = prominence_term * center_term
            ranked.append((score, lp, up, midpoint_error))
    ranked.sort(key=lambda row: row[0], reverse=True)
    score, lp, up, midpoint_error = ranked[0]
    alternatives = [
        {
            "lp_cm": row[1].wavenumber_cm,
            "up_cm": row[2].wavenumber_cm,
            "gap_cm": row[2].wavenumber_cm - row[1].wavenumber_cm,
            "midpoint_error_cm": row[3],
            "score": row[0],
        }
        for row in ranked[:5]
    ]
    return lp, up, score, alternatives


def fmt(value: float | None, digits: int = 2) -> str:
    return "not available" if value is None else f"{value:.{digits}f}"


def main() -> int:
    args = parse_args()
    omega_cm = args.omega_cm if args.omega_cm is not None else args.omega_au * AU_TO_CM
    analysis_range = tuple(sorted(args.analysis_range))
    normalization_range = (
        tuple(sorted(args.normalization_range))
        if args.normalization_range is not None
        else analysis_range
    )
    if args.smooth_fwhm_cm < 0:
        raise ValueError("--smooth-fwhm-cm cannot be negative")
    if args.min_rabi_gap_cm <= 0 or args.max_rabi_gap_cm <= args.min_rabi_gap_cm:
        raise ValueError(
            "Require 0 < --min-rabi-gap-cm < --max-rabi-gap-cm"
        )

    xc, yc, cavity_column = load_spectrum(args.cavity, args.intensity_column)
    xo, yo, outside_column = load_spectrum(args.outside, args.intensity_column)
    x, cavity_raw, outside_raw = make_common_grid(
        xc, yc, xo, yo, analysis_range
    )
    dx = float(np.median(np.diff(x)))
    cavity_smooth = smooth(cavity_raw, dx, args.smooth_fwhm_cm)
    outside_smooth = smooth(outside_raw, dx, args.smooth_fwhm_cm)
    cavity_norm, outside_norm, norm_scales = normalize(
        x,
        cavity_smooth,
        outside_smooth,
        normalization_range,
        args.normalization,
    )
    cavity_shape, outside_shape, _ = normalize(
        x, cavity_smooth, outside_smooth, normalization_range, "independent-max"
    )

    mask = window_mask(x, analysis_range)
    cavity_peaks = detect_peaks(
        x,
        cavity_smooth,
        mask,
        args.prominence_fraction,
        args.min_peak_distance_cm,
    )
    outside_peaks = detect_peaks(
        x,
        outside_smooth,
        mask,
        args.prominence_fraction,
        args.min_peak_distance_cm,
    )
    bare_peak = choose_bare_peak(
        outside_peaks, omega_cm, args.bare_search_halfwidth_cm
    )
    resolution_cm = (
        PS_TO_RESOLUTION_CM / args.trajectory_ps
        if args.trajectory_ps is not None and args.trajectory_ps > 0
        else None
    )
    exclusion_cm = (
        args.center_exclusion_cm
        if args.center_exclusion_cm is not None
        else max(args.smooth_fwhm_cm, resolution_cm or 0.0)
    )

    pair = None
    if bare_peak is not None:
        pair = choose_polariton_pair(
            cavity_peaks,
            omega_cm,
            bare_peak.wavenumber_cm,
            exclusion_cm,
            args.rabi_window_cm,
            args.midpoint_tolerance_cm,
            args.min_rabi_gap_cm,
            args.max_rabi_gap_cm,
        )

    result: dict[str, object] = {
        "input": {
            "cavity_file": str(args.cavity),
            "outside_file": str(args.outside),
            "cavity_intensity_column": cavity_column,
            "outside_intensity_column": outside_column,
        },
        "settings": {
            "omega_c_cm": omega_cm,
            "analysis_range_cm": analysis_range,
            "normalization_range_cm": normalization_range,
            "normalization": args.normalization,
            "normalization_scales": norm_scales,
            "smooth_fwhm_cm": args.smooth_fwhm_cm,
            "grid_spacing_cm": dx,
            "trajectory_ps": args.trajectory_ps,
            "fft_resolution_cm": resolution_cm,
            "center_exclusion_cm": exclusion_cm,
        },
        "bare_peak": asdict(bare_peak) if bare_peak else None,
        "cavity_peaks": [asdict(p) for p in cavity_peaks],
        "outside_peaks": [asdict(p) for p in outside_peaks],
        "assignment": None,
    }

    report: list[str] = []
    report.append("CAVITY / NO-CAVITY IR COMPARISON")
    report.append("=" * 39)
    report.append(f"Cavity frequency: {omega_cm:.3f} cm^-1")
    report.append(
        f"Smoothing: Gaussian FWHM = {args.smooth_fwhm_cm:.2f} cm^-1"
    )
    report.append(
        f"Normalization: {args.normalization}; scales = "
        f"{norm_scales['cavity']:.8g} (cavity), "
        f"{norm_scales['outside']:.8g} (outside)"
    )
    if resolution_cm is not None:
        report.append(
            f"Nominal physical FFT resolution from {args.trajectory_ps:g} ps: "
            f"{resolution_cm:.3f} cm^-1"
        )
    report.append("")

    assignment = None
    if bare_peak is None:
        report.append("No no-cavity peak could be assigned near the cavity frequency.")
        report.append("No Rabi gap was assigned.")
    elif pair is None:
        report.append(
            f"Assigned bare peak: {bare_peak.wavenumber_cm:.3f} cm^-1"
        )
        report.append("No cavity peak pair satisfying the LP/UP constraints was found.")
    else:
        lp, up, pair_score, alternatives = pair
        gap = up.wavenumber_cm - lp.wavenumber_cm
        detuning = omega_cm - bare_peak.wavenumber_cm
        resonant_gap = math.sqrt(max(gap * gap - detuning * detuning, 0.0))
        coupling_g = 0.5 * resonant_gap
        expected_midpoint = 0.5 * (omega_cm + bare_peak.wavenumber_cm)
        observed_midpoint = 0.5 * (lp.wavenumber_cm + up.wavenumber_cm)
        lp_ref, lp_ref_distance = nearest_peak(outside_peaks, lp.wavenumber_cm)
        up_ref, up_ref_distance = nearest_peak(outside_peaks, up.wavenumber_cm)
        center_mask = np.abs(x - expected_midpoint) <= exclusion_cm
        center_height = float(np.max(np.maximum(cavity_shape[center_mask], 0.0)))
        lp_i = int(np.argmin(np.abs(x - lp.wavenumber_cm)))
        up_i = int(np.argmin(np.abs(x - up.wavenumber_cm)))
        branch_height = 0.5 * (
            max(cavity_shape[lp_i], 0.0) + max(cavity_shape[up_i], 0.0)
        )
        center_to_branch = center_height / branch_height if branch_height > 0 else math.inf
        match_tolerance = max(args.smooth_fwhm_cm, resolution_cm or 0.0)
        warnings: list[str] = []
        if lp_ref_distance is not None and lp_ref_distance <= match_tolerance:
            warnings.append(
                f"LP lies {lp_ref_distance:.1f} cm^-1 from a no-cavity peak "
                f"at {lp_ref.wavenumber_cm:.1f} cm^-1"
            )
        if up_ref_distance is not None and up_ref_distance <= match_tolerance:
            warnings.append(
                f"UP lies {up_ref_distance:.1f} cm^-1 from a no-cavity peak "
                f"at {up_ref.wavenumber_cm:.1f} cm^-1"
            )
        if center_to_branch > 0.75:
            warnings.append(
                f"a strong central cavity feature remains "
                f"(central/mean-branch height = {center_to_branch:.2f})"
            )
        if resolution_cm is not None and gap < 3.0 * resolution_cm:
            warnings.append("the gap is less than three nominal FFT resolution elements")

        assignment = {
            "status": "provisional candidate; two IR spectra alone cannot prove hybridization",
            "lp_cm": lp.wavenumber_cm,
            "up_cm": up.wavenumber_cm,
            "observed_gap_cm": gap,
            "bare_vibration_cm": bare_peak.wavenumber_cm,
            "detuning_cm": detuning,
            "detuning_corrected_resonant_gap_cm": resonant_gap,
            "coupling_g_cm": coupling_g,
            "expected_midpoint_cm": expected_midpoint,
            "observed_midpoint_cm": observed_midpoint,
            "midpoint_error_cm": observed_midpoint - expected_midpoint,
            "pair_score": pair_score,
            "nearest_outside_peak_to_lp_cm": (
                lp_ref.wavenumber_cm if lp_ref else None
            ),
            "nearest_outside_peak_to_lp_distance_cm": lp_ref_distance,
            "nearest_outside_peak_to_up_cm": (
                up_ref.wavenumber_cm if up_ref else None
            ),
            "nearest_outside_peak_to_up_distance_cm": up_ref_distance,
            "central_to_mean_branch_height": center_to_branch,
            "warnings": warnings,
            "top_candidate_pairs": alternatives,
        }
        result["assignment"] = assignment

        report.extend(
            [
                f"Assigned bare peak: {bare_peak.wavenumber_cm:.3f} cm^-1",
                f"Candidate LP: {lp.wavenumber_cm:.3f} cm^-1",
                f"Candidate UP: {up.wavenumber_cm:.3f} cm^-1",
                f"Observed LP-UP gap: {gap:.3f} cm^-1",
                f"Detuning (omega_c - omega_v): {detuning:.3f} cm^-1",
                f"Detuning-corrected resonant gap: {resonant_gap:.3f} cm^-1",
                f"Inferred g = Omega_R/2: {coupling_g:.3f} cm^-1",
                f"Midpoint error: {observed_midpoint - expected_midpoint:.3f} cm^-1",
                "",
                "Assignment status: PROVISIONAL CANDIDATE.",
                "A matched q-spectrum or an avoided-crossing scan is required for confirmation.",
            ]
        )
        if warnings:
            report.append("")
            report.append("Cautions:")
            report.extend(f"- {warning}" for warning in warnings)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = prefix.with_suffix(".csv")
    dat_path = prefix.with_suffix(".dat")
    json_path = prefix.with_suffix(".json")
    txt_path = prefix.with_suffix(".txt")
    png_path = prefix.with_suffix(".png")
    table = np.column_stack(
        [
            x,
            cavity_raw,
            outside_raw,
            cavity_smooth,
            outside_smooth,
            cavity_norm,
            outside_norm,
            cavity_shape,
            outside_shape,
            cavity_norm - outside_norm,
        ]
    )
    header = (
        "wavenumber_cm-1,cavity_raw,outside_raw,cavity_smoothed,outside_smoothed,"
        "cavity_selected_normalized,outside_selected_normalized,cavity_independent_max,"
        "outside_independent_max,selected_normalized_difference"
    )
    np.savetxt(csv_path, table, delimiter=",", header=header, comments="")
    np.savetxt(dat_path, table, header=header.replace(",", " "), comments="# ")
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    txt_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    if not args.no_plot:
        try:
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(
                2, 1, figsize=(11, 8), sharex=True, layout="constrained"
            )
            ax = axes[0]
            ax.plot(x, outside_norm, color="#333333", lw=1.8, label="No cavity (lambda = 0)")
            ax.plot(x, cavity_norm, color="#1769aa", lw=1.8, label="Cavity")
            ax.axvline(omega_cm, color="#c62828", ls="--", lw=1.2, label="Cavity frequency")
            scale_label = (
                "shared reference scale"
                if args.normalization.startswith("reference-")
                else args.normalization.replace("-", " ")
            )
            ax.set_ylabel(f"IR intensity ({scale_label})")
            ax.legend(frameon=False, ncol=3, fontsize=9)
            ax.grid(alpha=0.18)

            ax2 = axes[1]
            ax2.plot(x, outside_shape, color="#333333", lw=1.7, label="No cavity")
            ax2.plot(x, cavity_shape, color="#1769aa", lw=1.7, label="Cavity")
            ax2.axvline(omega_cm, color="#c62828", ls="--", lw=1.2)
            if bare_peak is not None:
                ax2.axvline(bare_peak.wavenumber_cm, color="#777777", ls=":", lw=1.1)
            if assignment is not None:
                lp_cm = float(assignment["lp_cm"])
                up_cm = float(assignment["up_cm"])
                for value, label in ((lp_cm, "LP?"), (up_cm, "UP?")):
                    ax2.axvline(value, color="#ef6c00", ls="-.", lw=1.2)
                    ax2.text(
                        value,
                        1.02,
                        label,
                        color="#b45309",
                        ha="center",
                        va="bottom",
                        fontsize=9,
                    )
                ax2.set_title(
                    f"Candidate gap = {float(assignment['observed_gap_cm']):.1f} cm$^{{-1}}$ "
                    "(provisional)",
                    fontsize=11,
                )
            ax2.set_xlabel("Wavenumber / cm$^{-1}$")
            ax2.set_ylabel("Each spectrum / own maximum\n(shape comparison only)")
            ax2.set_xlim(*analysis_range)
            ax2.grid(alpha=0.18)
            fig.suptitle(
                f"Cavity versus no-cavity IR; Gaussian FWHM {args.smooth_fwhm_cm:g} cm$^{{-1}}$",
                fontsize=13,
            )
            fig.savefig(png_path, dpi=220)
            plt.close(fig)
        except ImportError:
            print("matplotlib is unavailable; skipped PNG output", file=sys.stderr)

    print("\n".join(report))
    print("")
    print(f"Wrote {dat_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {txt_path}")
    print(f"Wrote {json_path}")
    if not args.no_plot and png_path.exists():
        print(f"Wrote {png_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
