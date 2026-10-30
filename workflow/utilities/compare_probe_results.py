#!/usr/bin/env python3
"""Compare tight-reference and candidate eT probe JSON results."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics


def vector_difference(first: list[float], second: list[float]) -> list[float]:
    if len(first) != len(second):
        raise ValueError("Force arrays have different lengths")
    return [candidate - reference for reference, candidate in zip(first, second)]


def load(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    required = {"symbols", "energy_hartree", "forces_au", "metadata"}
    missing = required.difference(payload)
    if missing:
        raise ValueError(f"{path}: missing keys {sorted(missing)}")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="tight-reference probe JSON")
    parser.add_argument("candidate", help="candidate-profile probe JSON")
    parser.add_argument(
        "--max-force-rmse", type=float,
        help="optional acceptance limit for nuclear-force RMSE (Eh/bohr)",
    )
    parser.add_argument(
        "--max-force-error", type=float,
        help="optional acceptance limit for maximum nuclear-force error (Eh/bohr)",
    )
    parser.add_argument(
        "--max-polarized-dipole-error", type=float,
        help="optional acceptance limit for polarized-dipole error (ea0)",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    reference = load(Path(args.reference))
    candidate = load(Path(args.candidate))
    if reference["symbols"] != candidate["symbols"]:
        raise SystemExit("Reference and candidate atom symbols differ")

    symbols = list(reference["symbols"])
    if not symbols or symbols[-1] != "L":
        raise SystemExit("Expected the final probe particle to be L")
    nuclear_components = 3 * (len(symbols) - 1)
    differences = vector_difference(
        list(reference["forces_au"]), list(candidate["forces_au"])
    )[:nuclear_components]
    force_rmse = math.sqrt(statistics.fmean(value * value for value in differences))
    force_max = max(abs(value) for value in differences)

    reference_metadata = dict(reference["metadata"])
    candidate_metadata = dict(candidate["metadata"])
    reference_mu = float(reference_metadata["polarized_dipole_au"])
    candidate_mu = float(candidate_metadata["polarized_dipole_au"])
    dipole_error = abs(candidate_mu - reference_mu)
    energy_delta = float(candidate["energy_hartree"]) - float(
        reference["energy_hartree"]
    )
    reference_wall = float(reference_metadata["et_launcher_wall_seconds"])
    candidate_wall = float(candidate_metadata["et_launcher_wall_seconds"])

    print(f"Energy difference (candidate-reference): {energy_delta:+.10e} Eh")
    print(f"Nuclear-force RMSE: {force_rmse:.6e} Eh/bohr")
    print(f"Maximum nuclear-force error: {force_max:.6e} Eh/bohr")
    print(f"Polarized-dipole error: {dipole_error:.6e} ea0")
    print(f"Reference wall time: {reference_wall:.3f} s")
    print(f"Candidate wall time: {candidate_wall:.3f} s")
    print(f"Single-point speedup: {reference_wall / candidate_wall:.3f}x")

    failures = []
    if args.max_force_rmse is not None and force_rmse > args.max_force_rmse:
        failures.append("force RMSE")
    if args.max_force_error is not None and force_max > args.max_force_error:
        failures.append("maximum force error")
    if (
        args.max_polarized_dipole_error is not None
        and dipole_error > args.max_polarized_dipole_error
    ):
        failures.append("polarized-dipole error")
    if failures:
        print("REJECT: exceeded " + ", ".join(failures))
        return 2
    if any(
        value is not None
        for value in (
            args.max_force_rmse,
            args.max_force_error,
            args.max_polarized_dipole_error,
        )
    ):
        print("ACCEPT: all requested numerical limits passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
