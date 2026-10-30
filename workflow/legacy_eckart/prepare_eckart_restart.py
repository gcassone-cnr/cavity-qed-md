#!/usr/bin/env python3
"""Prepare an i-PI RESTART for an orientation-fixed H2O CBO-HF NVE run.

The default setup assumes the atom order O, H, H, L used by the eT/i-PI
CBO-HF calculations in this project.  It:

  1. rotates the H2O bisector onto the laboratory +z axis;
  2. puts the molecular plane in xz (H1-H2 projected along +x);
  3. applies the same rotation to the nuclear momenta;
  4. removes the H2O centre-of-mass momentum and rigid-body angular momentum;
  5. leaves photon q=Lx and p_q unchanged, while setting the unused Ly/Lz
     coordinates and momenta to zero; and
  6. replaces the i-PI motion block with Eckart-constrained NVE dynamics.

No external packages other than NumPy are required.  The original file is
never overwritten unless --overwrite is explicitly supplied.

Example
-------
    python3 prepare_eckart_restart.py RESTART \
        --output RESTART_eckart

Then launch i-PI with RESTART_eckart.  The Eckart reference is the aligned
geometry written by this script because no explicit reference values are put
in the constraint block.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np


FLOAT_RE = re.compile(
    r"[-+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[EeDd][-+]?\d+)?"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align H2O with the cavity polarization, remove its external "
            "momenta, and prepare an Eckart-constrained i-PI NVE restart."
        )
    )
    parser.add_argument("restart", type=Path, help="Input i-PI RESTART/XML file")
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output file (default: INPUTNAME_eckart)",
    )
    parser.add_argument(
        "--polarization", nargs=3, type=float, default=(0.0, 0.0, 1.0),
        metavar=("EX", "EY", "EZ"),
        help="Cavity-polarization vector (default: 0 0 1)",
    )
    parser.add_argument(
        "--plane-vector", nargs=3, type=float, default=(1.0, 0.0, 0.0),
        metavar=("PX", "PY", "PZ"),
        help=(
            "Direction for the projected H1-H2 vector; must not be parallel "
            "to the polarization (default: 1 0 0)"
        ),
    )
    parser.add_argument(
        "--water-indices", nargs=3, type=int, default=(0, 1, 2),
        metavar=("O", "H1", "H2"),
        help="Zero-based O,H,H atom indices (default: 0 1 2)",
    )
    parser.add_argument(
        "--photon-index", type=int, default=3,
        help="Zero-based photon pseudo-atom index (default: 3)",
    )
    parser.add_argument(
        "--timestep-fs", type=float,
        help="Override the timestep in the generated NVE block",
    )
    parser.add_argument(
        "--nsteps-geo", type=int, default=5,
        help="Geodesic substeps for constrained dynamics (default: 5)",
    )
    parser.add_argument(
        "--constraint-tolerance", type=float, default=1.0e-8,
        help="Constraint-solver tolerance (default: 1e-8)",
    )
    parser.add_argument(
        "--constraint-maxit", type=int, default=1000,
        help="Maximum constraint-solver iterations (default: 1000)",
    )
    parser.add_argument(
        "--no-motion-update", action="store_true",
        help="Transform q/p only; do not replace the motion block",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting the output file",
    )
    return parser.parse_args()


def element_match(text: str, tag: str) -> re.Match[str]:
    pattern = re.compile(
        rf"(?P<open><{re.escape(tag)}\b[^>]*>)(?P<body>.*?)(?P<close></{re.escape(tag)}>)",
        re.DOTALL,
    )
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"Could not find <{tag}>...</{tag}> in the restart")
    return match


def parse_numeric_element(text: str, tag: str) -> np.ndarray:
    body = element_match(text, tag).group("body")
    values = [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(body)]
    if not values:
        raise ValueError(f"No numerical values found in <{tag}>")
    return np.asarray(values, dtype=float)


def parse_names(text: str) -> list[str]:
    body = element_match(text, "names").group("body")
    body = body.replace("[", " ").replace("]", " ")
    return [item for item in re.split(r"[\s,]+", body.strip()) if item]


def replace_numeric_element(text: str, tag: str, values: np.ndarray) -> str:
    match = element_match(text, tag)
    flat = np.asarray(values, dtype=float).reshape(-1)
    chunks = []
    for start in range(0, flat.size, 4):
        chunks.append(", ".join(f"{x: .16e}" for x in flat[start : start + 4]))
    body = "\n          [ " + ",\n            ".join(chunks) + " ]\n         "
    return text[: match.start()] + match.group("open") + body + match.group("close") + text[match.end() :]


def unit(vector: np.ndarray, label: str) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1.0e-14:
        raise ValueError(f"Cannot normalize the near-zero {label}")
    return vector / norm


def skew(vector: np.ndarray) -> np.ndarray:
    x, y, z = vector
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def rotation_from_to(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Return a proper rotation that maps unit source onto unit target."""
    source = unit(source, "source direction")
    target = unit(target, "target direction")
    cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
    cross = np.cross(source, target)
    sine = float(np.linalg.norm(cross))

    if sine > 1.0e-14:
        k = skew(cross)
        return np.eye(3) + k + k @ k * ((1.0 - cosine) / (sine * sine))
    if cosine > 0.0:
        return np.eye(3)

    # Antiparallel vectors: choose any stable axis perpendicular to source.
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(source, trial))) > 0.8:
        trial = np.array([0.0, 1.0, 0.0])
    axis = unit(np.cross(source, trial), "180-degree rotation axis")
    return 2.0 * np.outer(axis, axis) - np.eye(3)


def rotation_about_axis(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = unit(axis, "rotation axis")
    k = skew(axis)
    return np.eye(3) + math.sin(angle) * k + (1.0 - math.cos(angle)) * (k @ k)


def bisector(positions: np.ndarray, indices: tuple[int, int, int]) -> np.ndarray:
    oxygen, hydrogen1, hydrogen2 = indices
    u1 = unit(positions[hydrogen1] - positions[oxygen], "O-H1 vector")
    u2 = unit(positions[hydrogen2] - positions[oxygen], "O-H2 vector")
    return unit(u1 + u2, "H-O-H bisector")


def alignment_rotation(
    positions: np.ndarray,
    indices: tuple[int, int, int],
    polarization: np.ndarray,
    plane_vector: np.ndarray,
) -> np.ndarray:
    polarization = unit(polarization, "polarization vector")
    plane_vector = plane_vector - np.dot(plane_vector, polarization) * polarization
    plane_vector = unit(plane_vector, "plane vector projected perpendicular to polarization")

    first = rotation_from_to(bisector(positions, indices), polarization)
    _, hydrogen1, hydrogen2 = indices
    hh = first @ (positions[hydrogen1] - positions[hydrogen2])
    hh -= np.dot(hh, polarization) * polarization
    hh = unit(hh, "projected H1-H2 vector")

    cosine = float(np.clip(np.dot(hh, plane_vector), -1.0, 1.0))
    sine = float(np.dot(polarization, np.cross(hh, plane_vector)))
    second = rotation_about_axis(polarization, math.atan2(sine, cosine))
    return second @ first


def water_com(positions: np.ndarray, masses: np.ndarray) -> np.ndarray:
    return np.sum(masses[:, None] * positions, axis=0) / np.sum(masses)


def total_angular_momentum(
    positions: np.ndarray, momenta: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    relative = positions - water_com(positions, masses)
    return np.sum(np.cross(relative, momenta), axis=0)


def remove_external_momenta(
    positions: np.ndarray, momenta: np.ndarray, masses: np.ndarray
) -> tuple[np.ndarray, float]:
    """Remove COM and rigid rotational momentum; return momenta and removed KE."""
    corrected = momenta.copy()
    kinetic_before = float(np.sum(corrected * corrected / masses[:, None]) / 2.0)

    total_mass = float(np.sum(masses))
    total_p = np.sum(corrected, axis=0)
    corrected -= masses[:, None] * total_p[None, :] / total_mass

    relative = positions - water_com(positions, masses)
    inertia = np.zeros((3, 3))
    identity = np.eye(3)
    for r_i, mass_i in zip(relative, masses):
        inertia += mass_i * (np.dot(r_i, r_i) * identity - np.outer(r_i, r_i))

    angular_momentum = np.sum(np.cross(relative, corrected), axis=0)
    angular_velocity = np.linalg.pinv(inertia, rcond=1.0e-13) @ angular_momentum
    corrected -= masses[:, None] * np.cross(angular_velocity, relative)

    # Remove any translation introduced at roundoff level.
    residual_p = np.sum(corrected, axis=0)
    corrected -= masses[:, None] * residual_p[None, :] / total_mass

    kinetic_after = float(np.sum(corrected * corrected / masses[:, None]) / 2.0)
    removed = kinetic_before - kinetic_after
    if abs(removed) < 1.0e-14:
        removed = 0.0
    return corrected, removed


def extract_timestep(text: str) -> tuple[float, str]:
    match = re.search(
        r"<timestep\b(?P<attrs>[^>]*)>(?P<value>.*?)</timestep>", text, re.DOTALL
    )
    if match is None:
        raise ValueError("Could not find a <timestep> in the input motion block")
    value_match = FLOAT_RE.search(match.group("value"))
    if value_match is None:
        raise ValueError("Could not parse the timestep value")
    units_match = re.search(r"units\s*=\s*['\"]([^'\"]+)['\"]", match.group("attrs"))
    units = units_match.group(1) if units_match else "atomic_unit"
    return float(value_match.group(0)), units


def replace_motion_block(
    text: str,
    water_indices: tuple[int, int, int],
    timestep: float,
    timestep_units: str,
    nsteps_geo: int,
    tolerance: float,
    maxit: int,
) -> str:
    pattern = re.compile(
        r"(?P<indent>^[ \t]*)<motion\b[^>]*>.*?</motion>", re.DOTALL | re.MULTILINE
    )
    match = pattern.search(text)
    if match is None:
        raise ValueError("Could not find the <motion>...</motion> block")
    indent = match.group("indent")
    i1 = indent + "   "
    i2 = indent + "      "
    i3 = indent + "         "
    atoms = ", ".join(str(i) for i in water_indices)
    replacement = (
        f"{indent}<motion mode='constrained_dynamics'>\n"
        f"{i1}<fixcom>False</fixcom>\n"
        f"{i1}<constrained_dynamics mode='nve'>\n"
        f"{i2}<timestep units='{timestep_units}'>{timestep:.16g}</timestep>\n"
        f"{i2}<nsteps_geo>{nsteps_geo}</nsteps_geo>\n"
        f"{i2}<csolver>\n"
        f"{i3}<maxit>{maxit}</maxit>\n"
        f"{i3}<norm_order>2</norm_order>\n"
        f"{i3}<tolerance>{tolerance:.8e}</tolerance>\n"
        f"{i2}</csolver>\n"
        f"{i2}<constraint mode='eckart'>\n"
        f"{i3}<atoms>[ {atoms} ]</atoms>\n"
        f"{i2}</constraint>\n"
        f"{i1}</constrained_dynamics>\n"
        f"{indent}</motion>"
    )
    return text[: match.start()] + replacement + text[match.end() :]


def main() -> int:
    args = parse_args()
    source = args.restart.resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Input restart does not exist: {source}")

    if args.output is None:
        if source.suffix:
            output = source.with_name(source.stem + "_eckart" + source.suffix)
        else:
            output = source.with_name(source.name + "_eckart")
    else:
        output = args.output.resolve()
    if output == source and not args.overwrite:
        raise FileExistsError("Refusing to overwrite the input; use a different --output")
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {output}; use --overwrite to replace it")

    text = source.read_text(encoding="utf-8")
    names = parse_names(text)
    natoms_match = re.search(r"<beads\b[^>]*\bnatoms=['\"](\d+)['\"]", text)
    nbeads_match = re.search(r"<beads\b[^>]*\bnbeads=['\"](\d+)['\"]", text)
    if natoms_match is None or nbeads_match is None:
        raise ValueError("Could not read natoms/nbeads from the <beads> tag")
    natoms = int(natoms_match.group(1))
    nbeads = int(nbeads_match.group(1))
    if nbeads != 1:
        raise ValueError("Eckart constrained dynamics is classical-only here; nbeads must be 1")
    if len(names) != natoms:
        raise ValueError(f"Found {len(names)} names but natoms={natoms}")

    q_flat = parse_numeric_element(text, "q")
    p_flat = parse_numeric_element(text, "p")
    masses = parse_numeric_element(text, "m")
    if q_flat.size != 3 * natoms or p_flat.size != 3 * natoms or masses.size != natoms:
        raise ValueError("The q, p, or m array size is inconsistent with natoms")
    positions = q_flat.reshape(natoms, 3)
    momenta = p_flat.reshape(natoms, 3)

    water_indices = tuple(args.water_indices)
    if len(set(water_indices)) != 3 or any(i < 0 or i >= natoms for i in water_indices):
        raise ValueError("--water-indices must contain three distinct valid atom indices")
    if args.photon_index < 0 or args.photon_index >= natoms:
        raise ValueError("--photon-index is outside the atom range")
    if args.photon_index in water_indices:
        raise ValueError("The photon index cannot also be a water-atom index")

    expected = ("O", "H", "H")
    actual = tuple(names[i].upper() for i in water_indices)
    if actual != expected:
        raise ValueError(
            f"Selected water atoms are {actual}, expected {expected}; check --water-indices"
        )

    polarization = unit(np.asarray(args.polarization, dtype=float), "polarization")
    plane_vector = np.asarray(args.plane_vector, dtype=float)
    water_pos = positions[list(water_indices)].copy()
    water_p = momenta[list(water_indices)].copy()
    water_m = masses[list(water_indices)].copy()

    initial_bisector = bisector(positions, water_indices)
    initial_total_p = np.sum(water_p, axis=0)
    initial_l = total_angular_momentum(water_pos, water_p, water_m)

    rotation = alignment_rotation(positions, water_indices, polarization, plane_vector)
    com = water_com(water_pos, water_m)
    water_pos = com + (rotation @ (water_pos - com).T).T
    water_p = (rotation @ water_p.T).T
    water_p, removed_ke = remove_external_momenta(water_pos, water_p, water_m)

    new_positions = positions.copy()
    new_momenta = momenta.copy()
    new_positions[list(water_indices)] = water_pos
    new_momenta[list(water_indices)] = water_p

    # q is Lx in this CBO-HF representation.  Preserve q and p_q exactly;
    # zero only the two unused pseudo-Cartesian directions.
    photon = args.photon_index
    q_photon_before = float(positions[photon, 0])
    p_photon_before = float(momenta[photon, 0])
    new_positions[photon, 1:] = 0.0
    new_momenta[photon, 1:] = 0.0

    final_bisector = bisector(new_positions, water_indices)
    final_total_p = np.sum(new_momenta[list(water_indices)], axis=0)
    final_l = total_angular_momentum(
        new_positions[list(water_indices)],
        new_momenta[list(water_indices)],
        water_m,
    )

    output_text = replace_numeric_element(text, "q", new_positions.reshape(-1))
    output_text = replace_numeric_element(output_text, "p", new_momenta.reshape(-1))

    if not args.no_motion_update:
        old_timestep, old_units = extract_timestep(text)
        if args.timestep_fs is None:
            timestep = old_timestep
            timestep_units = old_units
        else:
            timestep = args.timestep_fs
            timestep_units = "femtosecond"
        output_text = replace_motion_block(
            output_text,
            water_indices,
            timestep,
            timestep_units,
            args.nsteps_geo,
            args.constraint_tolerance,
            args.constraint_maxit,
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(output_text, encoding="utf-8")

    print(f"Input:  {source}")
    print(f"Output: {output}")
    print(f"Atoms:  {', '.join(f'{i}:{names[i]}' for i in range(natoms))}")
    print(f"Initial bisector dot polarization: {np.dot(initial_bisector, polarization): .12e}")
    print(f"Final   bisector dot polarization: {np.dot(final_bisector, polarization): .12e}")
    print(f"Initial H2O total momentum norm:    {np.linalg.norm(initial_total_p): .12e} a.u.")
    print(f"Final   H2O total momentum norm:    {np.linalg.norm(final_total_p): .12e} a.u.")
    print(f"Initial H2O angular momentum norm:  {np.linalg.norm(initial_l): .12e} a.u.")
    print(f"Final   H2O angular momentum norm:  {np.linalg.norm(final_l): .12e} a.u.")
    print(f"Removed external kinetic energy:   {removed_ke: .12e} hartree")
    print(f"Photon q=Lx preserved:              {q_photon_before: .16e}")
    print(f"Photon p_q preserved:               {p_photon_before: .16e}")
    print("Unused Ly/Lz positions and momenta: set to exactly zero")
    if args.no_motion_update:
        print("Motion block: unchanged (--no-motion-update)")
    else:
        print("Motion block: Eckart-constrained NVE; fixed Cartesian DOFs removed")
    print("\nInspect the checks above before launching i-PI.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
