#!/usr/bin/env python3
"""eT fixed-displacement QED-HF client for the base i-PI socket protocol.

The i-PI system contains the nuclei followed by one pseudo-atom ``L``.  The
x coordinate of L is the mass-weighted cavity coordinate q in atomic units;
its y and z coordinates must remain zero. Ordinary MD freezes these DOFs;
Eckart NVE starts them with zero momentum and the driver returns zero force.

At each geometry this client converts q to eT's coherent displacement
z = q*sqrt(omega/2), runs fixed-displacement QED-HF, returns the nuclear
forces, and evaluates the photon force

    F_q = -omega**2*q - omega*lambda*<mu.e>.

The patched analytic eT gradient is the normal mode. A slow central
finite-difference nuclear-gradient fallback is available for validation.
Each client can reuse its last converged SCF orbitals with an automatic SAD
fallback.  A rolling calculation directory prevents permanent per-step output
growth while preserving failed attempts and optional sparse snapshots.
Only the unbatched base i-PI protocol is implemented.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import struct
import subprocess
import sys
import tempfile
import time
from typing import Iterable, Sequence


HEADER_LENGTH = 12
INT32 = struct.Struct("=i")
FLOAT64 = struct.Struct("=d")
BOHR_PER_ANGSTROM = 1.8897261254578281
DEBYE_TO_AU = 0.393430270
FLOAT_TOKEN = (
    r"[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))"
    r"(?:[EeDd][+-]?\d+)?"
)
AXES = {"x": 0, "y": 1, "z": 2}


class PeerClosed(ConnectionError):
    """The peer closed between protocol messages."""


def header(message: str) -> bytes:
    return message.upper().ljust(HEADER_LENGTH).encode("ascii")


def recv_exact(sock: socket.socket, count: int) -> bytes:
    chunks = bytearray()
    while len(chunks) < count:
        chunk = sock.recv(count - len(chunks))
        if not chunk:
            if not chunks:
                raise PeerClosed("i-PI closed the socket")
            raise ConnectionError("i-PI closed the socket during a protocol message")
        chunks.extend(chunk)
    return bytes(chunks)


def recv_header(sock: socket.socket) -> str:
    return recv_exact(sock, HEADER_LENGTH).decode("ascii").strip().upper()


def recv_int32(sock: socket.socket) -> int:
    return INT32.unpack(recv_exact(sock, INT32.size))[0]


def recv_float64s(sock: socket.socket, count: int) -> tuple[float, ...]:
    data = recv_exact(sock, FLOAT64.size * count)
    return struct.unpack(f"={count}d", data)


def parse_number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def parse_xyz(path: Path) -> tuple[list[str], list[tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"{path} is empty")
    natoms = int(lines[0].strip())
    atom_lines = lines[2 : 2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError(f"{path} declares {natoms} atoms but has fewer coordinates")
    symbols: list[str] = []
    positions: list[tuple[float, float, float]] = []
    for line in atom_lines:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Invalid XYZ atom line: {line!r}")
        symbols.append(fields[0])
        positions.append(tuple(float(value) for value in fields[1:4]))
    return symbols, positions


def flatten(rows: Iterable[Sequence[float]]) -> list[float]:
    return [value for row in rows for value in row]


def normalize_polarization(values: Sequence[float]) -> tuple[float, float, float]:
    if len(values) != 3:
        raise ValueError("Polarization needs exactly three components")
    norm = math.sqrt(sum(value * value for value in values))
    if norm == 0.0:
        raise ValueError("Polarization vector cannot be zero")
    return tuple(value / norm for value in values)  # type: ignore[return-value]


def render_input(
    template: str,
    symbols: Sequence[str],
    positions_bohr: Sequence[float],
    frequency: float,
    coupling: float,
    polarization: Sequence[float],
    q: float,
    include_gradient: bool = True,
    scf_restart: bool = False,
) -> str:
    if len(positions_bohr) != 3 * len(symbols):
        raise ValueError("Position array and nuclear symbol list have inconsistent sizes")
    geometry = []
    for atom, symbol in enumerate(symbols):
        x, y, z = positions_bohr[3 * atom : 3 * atom + 3]
        geometry.append(f"  {symbol:<3s} {x: .14f} {y: .14f} {z: .14f}")
    coherent_z = q * math.sqrt(0.5 * frequency)
    replacements = {
        "@@FREQUENCY@@": f"{frequency:.16g}",
        "@@COUPLING@@": f"{coupling:.16g}",
        "@@POL_X@@": f"{polarization[0]:.16g}",
        "@@POL_Y@@": f"{polarization[1]:.16g}",
        "@@POL_Z@@": f"{polarization[2]:.16g}",
        "@@COHERENT_Z@@": f"{coherent_z:.16g}",
        "@@GEOMETRY_BOHR@@": "\n".join(geometry),
    }
    rendered = template
    for marker, value in replacements.items():
        rendered = rendered.replace(marker, value)
    if not include_gradient:
        rendered = "\n".join(
            line
            for line in rendered.splitlines()
            if line.strip().lower() != "molecular gradient"
        ) + ("\n" if rendered.endswith("\n") else "")
    if scf_restart:
        rendered = enable_scf_restart(rendered)
    leftovers = sorted(set(re.findall(r"@@[A-Z0-9_]+@@", rendered)))
    if leftovers:
        raise ValueError(f"Unreplaced template markers: {', '.join(leftovers)}")
    return rendered


def enable_scf_restart(text: str) -> str:
    """Insert eT's SCF restart keyword without changing other sections."""
    lines = text.splitlines()
    section_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.strip().lower() == "- solver scf"
        ),
        None,
    )
    if section_index is None:
        raise ValueError("Cannot enable SCF restart: template has no '- solver scf'")
    section_end = next(
        (
            index
            for index in range(section_index + 1, len(lines))
            if lines[index].lstrip().startswith("-")
        ),
        len(lines),
    )
    if any(
        line.strip().lower() == "restart"
        for line in lines[section_index + 1 : section_end]
    ):
        return text
    lines.insert(section_index + 1, "  restart")
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


def parse_scf_iterations(text: str) -> int | None:
    matches = re.findall(
        r"Convergence\s+criterion\s+met\s+in\s+(\d+)\s+iterations",
        text,
        flags=re.IGNORECASE,
    )
    return int(matches[-1]) if matches else None


def parse_energy(text: str, custom_regex: str | None = None) -> float:
    patterns = []
    if custom_regex:
        patterns.append(re.compile(custom_regex, re.IGNORECASE | re.MULTILINE))
    patterns.extend(
        [
            re.compile(
                rf"(?im)^[^\n]*(?:QED[\s-]*HF)[^\n]*?\benergy\b"
                rf"[^\d+\-.]*({FLOAT_TOKEN})"
            ),
            re.compile(
                rf"(?im)^\s*(?:total\s+)?(?:Hartree[\s-]*Fock|HF)"
                rf"[^\n]*?\benergy\b[^\d+\-.]*({FLOAT_TOKEN})"
            ),
            re.compile(
                rf"(?im)^\s*total\s+(?:electronic\s+)?energy"
                rf"[^\d+\-.]*({FLOAT_TOKEN})"
            ),
        ]
    )
    for pattern in patterns:
        matches = list(pattern.finditer(text))
        if matches:
            return parse_number(matches[-1].group(1))
    raise ValueError(
        "Could not identify the final QED-HF energy in the eT main output. "
        "Pass --energy-regex with one capturing group if this eT version "
        "uses a different label."
    )


def parse_total_dipole_debye(text: str) -> tuple[float, float, float]:
    marker = text.rfind("Dipole moment [Debye]")
    if marker < 0:
        raise ValueError("Could not find eT's 'Dipole moment [Debye]' table")
    component_pattern = re.compile(
        rf"(?im)^\s*([xyz])\s+({FLOAT_TOKEN})\s+({FLOAT_TOKEN})\s+({FLOAT_TOKEN})\s*$"
    )
    components: dict[str, float] = {}
    for match in component_pattern.finditer(text[marker:]):
        components[match.group(1).lower()] = parse_number(match.group(4))
        if len(components) == 3:
            break
    if set(components) != set(AXES):
        raise ValueError("Could not parse all total-dipole components from eT output")
    return components["x"], components["y"], components["z"]


def numeric_rows(text: str) -> list[tuple[float, float, float]]:
    """Return the last three numbers from lines containing at least three."""
    rows: list[tuple[float, float, float]] = []
    for line in text.splitlines():
        tokens = re.findall(FLOAT_TOKEN, line)
        if len(tokens) >= 3:
            values = tuple(parse_number(token) for token in tokens[-3:])
            rows.append(values)  # type: ignore[arg-type]
    return rows


def parse_gradient_file(path: Path, natoms: int) -> list[float]:
    rows = numeric_rows(path.read_text(encoding="utf-8", errors="replace"))
    if len(rows) < natoms:
        raise ValueError(f"{path.name} contains fewer than {natoms} numeric rows")
    gradient = flatten(rows[-natoms:])
    if len(gradient) != 3 * natoms or not all(math.isfinite(x) for x in gradient):
        raise ValueError(f"Invalid gradient parsed from {path.name}")
    if max(abs(x) for x in gradient) > 10.0:
        raise ValueError(
            f"Implausibly large gradient parsed from {path.name}; "
            "check its column layout"
        )
    return gradient


def find_gradient(
    run_directory: Path, natoms: int, gradient_glob: str
) -> tuple[list[float], Path]:
    candidates = sorted(
        (
            path
            for path in run_directory.rglob(gradient_glob)
            if path.is_file()
            and path.suffix.lower() not in {".inp", ".xml", ".xyz"}
        ),
        key=lambda path: path.stat().st_mtime_ns,
        reverse=True,
    )
    failures = []
    for candidate in candidates:
        try:
            return parse_gradient_file(candidate, natoms), candidate
        except ValueError as exc:
            failures.append(str(exc))
    names = ", ".join(path.name for path in candidates) or "none"
    details = "; ".join(failures[:3])
    raise ValueError(
        f"Could not parse a {natoms}-atom gradient. Candidate files: {names}. "
        f"{details} Set --gradient-glob for this eT version if needed."
    )


@dataclass
class ETPoint:
    energy: float
    dipole_au: tuple[float, float, float]
    gradient: list[float] | None
    gradient_file: str | None
    scf_iterations: int | None = None
    launcher_wall_seconds: float | None = None
    restart_loaded: bool = False
    restart_fallback: bool = False


class ETEvaluator:
    def __init__(self, args: argparse.Namespace, nuclear_symbols: Sequence[str]):
        self.et_launch = args.et_launch
        self.et_binary = getattr(args, "et_binary", None)
        self.basis_dir = getattr(args, "basis_dir", None)
        if self.et_binary:
            self.et_binary = str(Path(self.et_binary).resolve())
            if not os.access(self.et_binary, os.X_OK):
                raise ValueError(f"Cannot execute {self.et_binary}")
            if not self.basis_dir or not Path(self.basis_dir).is_dir():
                raise ValueError("Direct execution requires --basis-dir")
        self.template_path = Path(args.template).resolve()
        self.template = self.template_path.read_text(encoding="utf-8")
        method = re.search(r"(?ims)^\s*-\s*method\s*\n(.*?)(?=^\s*-\s*\w|\Z)", self.template)
        if not method or method.group(1).strip().lower() != "qed-hf":
            raise ValueError("This CBO driver requires method qed-hf; CBO-DFT forces are not implemented")
        self.workdir = Path(args.workdir).resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.nuclear_symbols = list(nuclear_symbols)
        self.frequency = args.frequency
        self.coupling = args.coupling
        self.polarization = normalize_polarization(args.polarization)
        self.photon_axis = AXES[args.photon_axis]
        self.threads = args.threads
        self.energy_regex = args.energy_regex
        self.gradient_glob = args.gradient_glob
        self.gradient_mode = args.gradient_mode
        self.fd_nuclear_step = args.fd_nuclear_step
        self.scf_restart = getattr(args, "scf_restart", "auto")
        self.keep_every = getattr(args, "keep_every", 0)
        self.failure_snapshots = not getattr(args, "no_failure_snapshots", False)
        self.current_directory = self.workdir / "current"
        self.restart_directory = self.workdir / "scf_restart"
        self.restart_stage = self.workdir / "scf_restart_stage"
        self.snapshot_directory = self.workdir / "snapshots"
        self.failure_directory = self.workdir / "failures"
        self.call_index = 0
        self.step = 0

    @property
    def natoms(self) -> int:
        return len(self.nuclear_symbols) + 1

    def _restart_available(self) -> bool:
        return self.restart_directory.is_dir() and any(
            self.restart_directory.glob("orbital_*")
        )

    @staticmethod
    def _replace_directory(path: Path) -> None:
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)

    def _archive_failure(self, label: str, kind: str) -> Path | None:
        if not self.failure_snapshots:
            return None
        safe_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)
        destination = (
            self.failure_directory
            / f"call_{self.call_index:06d}_{safe_label}_{kind}"
        )
        serial = 1
        while destination.exists():
            destination = (
                self.failure_directory
                / f"call_{self.call_index:06d}_{safe_label}_{kind}_{serial}"
            )
            serial += 1
        destination.mkdir(parents=True)
        if self.current_directory.exists():
            shutil.copytree(
                self.current_directory,
                destination / "calculation",
                dirs_exist_ok=True,
            )
        if self.restart_stage.exists():
            shutil.copytree(
                self.restart_stage,
                destination / "restart_stage",
                dirs_exist_ok=True,
            )
        return destination

    def _promote_restart(self) -> None:
        saved = self.restart_stage / "et_restart_files"
        if not saved.is_dir() or not any(saved.glob("orbital_*")):
            raise RuntimeError(
                "eT completed but did not save orbital restart files in "
                f"{saved}"
            )
        backup = self.workdir / "scf_restart_previous"
        if backup.exists():
            shutil.rmtree(backup)
        had_previous = self.restart_directory.exists()
        if had_previous:
            os.replace(self.restart_directory, backup)
        try:
            os.replace(saved, self.restart_directory)
        except Exception:
            if had_previous and backup.exists():
                os.replace(backup, self.restart_directory)
            raise
        if backup.exists():
            shutil.rmtree(backup)
        if self.restart_stage.exists():
            shutil.rmtree(self.restart_stage)

    def _write_evaluation_record(
        self,
        label: str,
        point: ETPoint,
        q: float,
        need_gradient: bool,
    ) -> None:
        record = {
            "call_index": self.call_index,
            "label": label,
            "q_au": q,
            "need_gradient": need_gradient,
            "energy_hartree": point.energy,
            "scf_iterations": point.scf_iterations,
            "launcher_wall_seconds": point.launcher_wall_seconds,
            "restart_loaded": point.restart_loaded,
            "restart_fallback": point.restart_fallback,
        }
        (self.current_directory / "evaluation.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def _run_et_once(
        self,
        label: str,
        nuclear_positions_bohr: Sequence[float],
        q: float,
        need_gradient: bool,
        load_restart: bool,
    ) -> ETPoint:
        self._replace_directory(self.current_directory)
        if self.scf_restart == "auto":
            self._replace_directory(self.restart_stage)
        input_path = self.current_directory / ("eT.inp" if self.et_binary else "et.inp")
        output_path = self.current_directory / ("eT.out" if self.et_binary else "et.out")
        input_path.write_text(
            render_input(
                self.template,
                self.nuclear_symbols,
                nuclear_positions_bohr,
                self.frequency,
                self.coupling,
                self.polarization,
                q,
                include_gradient=need_gradient,
                scf_restart=load_restart,
            ),
            encoding="utf-8",
        )
        command = [
            self.et_launch,
            str(input_path),
            "--omp",
            str(self.threads),
            "-of",
            str(output_path),
            "-nt",
            "-s",
        ]
        if self.scf_restart == "auto":
            command.extend(["-save", str(self.restart_stage)])
            if load_restart:
                command.extend(["-load", str(self.restart_directory)])
        environment = os.environ.copy()
        if self.et_binary:
            command = [self.et_binary]
            environment["AO_BASIS_PATH"] = str(Path(self.basis_dir).resolve())
            environment["LIBINT_DATA_PATH"] = environment["AO_BASIS_PATH"]
            environment["OMP_NUM_THREADS"] = str(self.threads)
            if load_restart:
                for orbital in self.restart_directory.glob("orbital_*"):
                    shutil.copy2(orbital, self.current_directory / orbital.name)
        started = time.perf_counter()
        completed = subprocess.run(
            command,
            cwd=self.current_directory,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )
        launcher_wall_seconds = time.perf_counter() - started
        (self.current_directory / "launcher.stdout").write_text(
            completed.stdout, encoding="utf-8"
        )
        (self.current_directory / "launcher.stderr").write_text(
            completed.stderr, encoding="utf-8"
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"eT failed with exit status {completed.returncode}; "
                f"inspect {self.current_directory}"
            )
        if not output_path.exists():
            raise FileNotFoundError(f"eT did not create {output_path}")
        output_text = output_path.read_text(encoding="utf-8", errors="replace")
        if "eT terminated successfully!" not in output_text:
            raise RuntimeError(
                f"eT did not report successful termination; inspect {output_path}"
            )
        energy = parse_energy(output_text, self.energy_regex)
        dipole_debye = parse_total_dipole_debye(output_text)
        dipole_au = tuple(value * DEBYE_TO_AU for value in dipole_debye)
        if not all(math.isfinite(v) for v in (energy, *dipole_au)):
            raise ValueError("Nonfinite electronic energy or dipole")
        gradient = None
        gradient_file = None
        if need_gradient:
            gradient, gradient_path = find_gradient(
                self.current_directory,
                len(self.nuclear_symbols),
                self.gradient_glob,
            )
            gradient_file = gradient_path.name
        point = ETPoint(
            energy=energy,
            dipole_au=dipole_au,
            gradient=gradient,
            gradient_file=gradient_file,
            scf_iterations=parse_scf_iterations(output_text),
            launcher_wall_seconds=launcher_wall_seconds,
            restart_loaded=load_restart,
        )
        self._write_evaluation_record(label, point, q, need_gradient)
        if self.scf_restart == "auto":
            if self.et_binary:
                saved = self.restart_stage / "et_restart_files"
                saved.mkdir()
                for orbital in self.current_directory.glob("orbital_*"):
                    shutil.copy2(orbital, saved / orbital.name)
            self._promote_restart()
        return point

    def _run_et(
        self,
        calculation_directory: Path,
        nuclear_positions_bohr: Sequence[float],
        q: float,
        need_gradient: bool,
    ) -> ETPoint:
        label = calculation_directory.name
        load_restart = self.scf_restart == "auto" and self._restart_available()
        try:
            try:
                return self._run_et_once(
                    label,
                    nuclear_positions_bohr,
                    q,
                    need_gradient,
                    load_restart,
                )
            except Exception as restart_error:
                failure_kind = "restart_failed" if load_restart else "cold_failed"
                archived = self._archive_failure(label, failure_kind)
                if not load_restart:
                    raise
                location = f"; snapshot: {archived}" if archived else ""
                print(
                    "WARNING: restarted SCF failed; retrying this geometry from SAD"
                    f"{location}. Cause: {restart_error}",
                    file=sys.stderr,
                    flush=True,
                )
                try:
                    point = self._run_et_once(
                        label,
                        nuclear_positions_bohr,
                        q,
                        need_gradient,
                        False,
                    )
                    point.restart_fallback = True
                    self._write_evaluation_record(label, point, q, need_gradient)
                    return point
                except Exception:
                    self._archive_failure(label, "sad_fallback_failed")
                    raise
        finally:
            self.call_index += 1

    def _finite_difference_gradient(
        self,
        step_directory: Path,
        nuclear_positions_bohr: Sequence[float],
        q: float,
    ) -> list[float]:
        gradient = []
        for coordinate in range(len(nuclear_positions_bohr)):
            plus = list(nuclear_positions_bohr)
            minus = list(nuclear_positions_bohr)
            plus[coordinate] += self.fd_nuclear_step
            minus[coordinate] -= self.fd_nuclear_step
            plus_point = self._run_et(
                step_directory / f"fd_{coordinate:03d}_plus", plus, q, False
            )
            minus_point = self._run_et(
                step_directory / f"fd_{coordinate:03d}_minus", minus, q, False
            )
            gradient.append(
                (plus_point.energy - minus_point.energy)
                / (2.0 * self.fd_nuclear_step)
            )
        return gradient

    def _snapshot_current(self, completed_count: int) -> Path | None:
        if (
            self.keep_every <= 0
            or completed_count % self.keep_every != 0
            or not self.current_directory.exists()
        ):
            return None
        destination = self.snapshot_directory / f"force_{completed_count:06d}"
        if destination.exists():
            raise FileExistsError(
                f"Sparse snapshot already exists: {destination}; "
                "use a new --workdir or remove the conflicting snapshot"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(self.current_directory, destination)
        return destination

    def evaluate_energy(self, positions_bohr: Sequence[float]) -> float:
        """Evaluate energy only, used by force checks and FD validation."""
        if len(positions_bohr) != 3 * self.natoms:
            raise ValueError(
                f"Received {len(positions_bohr)} coordinates; "
                f"expected {3 * self.natoms}"
            )
        nuclear_size = 3 * len(self.nuclear_symbols)
        nuclear_positions = list(positions_bohr[:nuclear_size])
        photon_position = positions_bohr[nuclear_size : nuclear_size + 3]
        q = photon_position[self.photon_axis]
        return self._run_et(
            Path(f"force_check_{self.call_index:06d}"),
            nuclear_positions,
            q,
            False,
        ).energy

    def evaluate(
        self, positions_bohr: Sequence[float]
    ) -> tuple[float, list[float], dict[str, object]]:
        if len(positions_bohr) != 3 * self.natoms:
            raise ValueError(
                f"Received {len(positions_bohr)} coordinates; "
                f"expected {3 * self.natoms}"
            )
        step_directory = Path(f"step_{self.step:06d}")
        nuclear_size = 3 * len(self.nuclear_symbols)
        nuclear_positions = list(positions_bohr[:nuclear_size])
        photon_position = positions_bohr[nuclear_size : nuclear_size + 3]
        q = photon_position[self.photon_axis]

        reference = self._run_et(
            step_directory / "reference",
            nuclear_positions,
            q,
            self.gradient_mode == "analytic",
        )
        if self.gradient_mode == "analytic":
            if reference.gradient is None:
                raise RuntimeError("Internal error: analytic gradient was not parsed")
            nuclear_gradient = reference.gradient
        else:
            nuclear_gradient = self._finite_difference_gradient(
                step_directory, nuclear_positions, q
            )

        polarized_dipole = sum(
            component * direction
            for component, direction in zip(reference.dipole_au, self.polarization)
        )
        photon_force = -(
            self.frequency**2 * q
            + self.frequency * self.coupling * polarized_dipole
        )
        forces = [-value for value in nuclear_gradient] + [0.0, 0.0, 0.0]
        forces[nuclear_size + self.photon_axis] = photon_force
        metadata: dict[str, object] = {
            "execution_backend": "direct" if self.et_binary else "launcher",
            "et_step": self.step,
            "energy_hartree": reference.energy,
            "q_au": q,
            "coherent_z": q * math.sqrt(0.5 * self.frequency),
            "dipole_au": reference.dipole_au,
            "polarized_dipole_au": polarized_dipole,
            "photon_force_au": photon_force,
            "nuclear_gradient_mode": self.gradient_mode,
            "gradient_file": reference.gradient_file,
            "scf_iterations": reference.scf_iterations,
            "et_launcher_wall_seconds": reference.launcher_wall_seconds,
            "scf_restart_loaded": reference.restart_loaded,
            "scf_restart_fallback": reference.restart_fallback,
        }
        snapshot = self._snapshot_current(self.step + 1)
        metadata["calculation_directory"] = str(self.current_directory)
        metadata["sparse_snapshot"] = str(snapshot) if snapshot else None
        self.step += 1
        return reference.energy, forces, metadata


def connect_unix(name: str, timeout: float) -> socket.socket:
    path = name if name.startswith("/") else f"/tmp/ipi_{name}"
    deadline = time.monotonic() + timeout
    last_error: OSError | None = None
    while time.monotonic() < deadline:
        sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        try:
            sock.connect(path)
            return sock
        except OSError as exc:
            last_error = exc
            sock.close()
            time.sleep(0.2)
    raise TimeoutError(f"Could not connect to i-PI UNIX socket {path}: {last_error}")


def connect_inet(host: str, port: int, timeout: float) -> socket.socket:
    deadline = time.monotonic() + timeout
    last_error: OSError | None = None
    while time.monotonic() < deadline:
        try:
            return socket.create_connection((host, port), timeout=5.0)
        except OSError as exc:
            last_error = exc
            time.sleep(0.2)
    raise TimeoutError(f"Could not connect to i-PI at {host}:{port}: {last_error}")


def send_force(
    sock: socket.socket,
    energy: float,
    forces: Sequence[float],
    metadata: dict[str, object],
) -> None:
    natoms = len(forces) // 3
    extra = json.dumps(metadata, separators=(",", ":")).encode("utf-8")
    payload = bytearray(header("FORCEREADY"))
    payload.extend(FLOAT64.pack(energy))
    payload.extend(INT32.pack(natoms))
    payload.extend(struct.pack(f"={3 * natoms}d", *forces))
    payload.extend(struct.pack("=9d", *([0.0] * 9)))
    payload.extend(INT32.pack(len(extra)))
    payload.extend(extra)
    sock.sendall(payload)


def serve(sock: socket.socket, evaluator: ETEvaluator) -> None:
    have_data = False
    result: tuple[float, list[float], dict[str, object]] | None = None
    print("Connected to i-PI; waiting for coordinates", flush=True)
    while True:
        try:
            command = recv_header(sock)
        except PeerClosed:
            if have_data:
                raise ConnectionError("i-PI closed before collecting the pending force")
            # i-PI 3.2.0 closes the connection without an EXIT message.
            # launch.py separately verifies the server's successful exit.
            print("i-PI closed the idle socket", flush=True)
            return
        if command == "STATUS":
            sock.sendall(header("HAVEDATA" if have_data else "READY"))
        elif command == "POSDATA":
            recv_float64s(sock, 9)
            recv_float64s(sock, 9)
            natoms = recv_int32(sock)
            if natoms != evaluator.natoms:
                raise ValueError(
                    f"i-PI sent {natoms} particles; expected {evaluator.natoms}"
                )
            positions = recv_float64s(sock, 3 * natoms)
            result = evaluator.evaluate(positions)
            have_data = True
            print(
                f"force evaluation {evaluator.step - 1}: "
                f"E = {result[0]:.14f} Eh, q = "
                f"{result[2]['q_au']:.8e} au, "
                f"SCF iterations = {result[2]['scf_iterations']}, "
                f"restart = {result[2]['scf_restart_loaded']}, "
                f"wall = {result[2]['et_launcher_wall_seconds']:.3f} s",
                flush=True,
            )
        elif command == "GETFORCE":
            if not have_data or result is None:
                raise RuntimeError("i-PI requested forces before an evaluation")
            send_force(sock, *result)
            have_data = False
            result = None
        elif command == "EXIT":
            print("i-PI requested exit", flush=True)
            return
        else:
            raise RuntimeError(f"Unexpected i-PI command {command!r}")


def parse_force_check(value: str, natoms: int) -> tuple[int, int]:
    atom_text, axis_text = value.split(",", maxsplit=1)
    atom = int(atom_text)
    axis_key = axis_text.lower()
    axis = AXES[axis_key] if axis_key in AXES else int(axis_text)
    if not (0 <= atom < natoms and 0 <= axis <= 2):
        raise ValueError("--check-force must look like ATOM,AXIS, e.g. 0,y")
    return atom, axis


def run_probe(args: argparse.Namespace, evaluator: ETEvaluator) -> int:
    symbols, positions_angstrom = parse_xyz(Path(args.probe))
    expected_symbols = evaluator.nuclear_symbols + ["L"]
    if symbols != expected_symbols:
        raise ValueError(
            f"Probe XYZ symbols {symbols} disagree with expected {expected_symbols}"
        )
    # XYZ stores nuclei in angstrom and the final L.x as q in atomic units.
    positions_bohr = [
        value * BOHR_PER_ANGSTROM for value in flatten(positions_angstrom[:-1])
    ] + list(positions_angstrom[-1])
    if args.probe_q is not None:
        photon_offset = 3 * len(evaluator.nuclear_symbols)
        positions_bohr[photon_offset + evaluator.photon_axis] = args.probe_q
    energy, forces, metadata = evaluator.evaluate(positions_bohr)
    print(f"Probe energy: {energy:.14f} Eh")
    print(f"Nuclear gradient mode: {metadata['nuclear_gradient_mode']}")
    print(f"Photon coordinate q: {metadata['q_au']:.10e} au")
    print(f"Polarized dipole: {metadata['polarized_dipole_au']:.10e} ea0")
    print("Forces (atomic units):")
    for symbol, vector in zip(symbols, zip(*[iter(forces)] * 3)):
        print(
            f"  {symbol:<3s} {vector[0]: .10e} "
            f"{vector[1]: .10e} {vector[2]: .10e}"
        )

    failed = False
    force_checks: list[dict[str, float | int]] = []
    for check in args.check_force:
        atom, axis = parse_force_check(check, len(symbols))
        coordinate = 3 * atom + axis
        plus = positions_bohr.copy()
        minus = positions_bohr.copy()
        plus[coordinate] += args.fd_check_step
        minus[coordinate] -= args.fd_check_step
        plus_energy = evaluator.evaluate_energy(plus)
        minus_energy = evaluator.evaluate_energy(minus)
        fd_force = -(plus_energy - minus_energy) / (2.0 * args.fd_check_step)
        reported_force = forces[coordinate]
        difference = abs(fd_force - reported_force)
        print(
            f"Force check atom {atom} axis {axis}: "
            f"reported={reported_force:.10e}, "
            f"finite-difference={fd_force:.10e}, "
            f"|delta|={difference:.3e} Eh/bohr"
        )
        force_checks.append(
            {
                "atom": atom,
                "axis": axis,
                "reported_force_au": reported_force,
                "finite_difference_force_au": fd_force,
                "absolute_difference_au": difference,
            }
        )
        if difference > args.fd_tolerance:
            failed = True
    if args.probe_json:
        probe_record = {
            "probe": str(Path(args.probe).resolve()),
            "template": str(Path(args.template).resolve()),
            "threads": args.threads,
            "symbols": symbols,
            "energy_hartree": energy,
            "forces_au": forces,
            "metadata": metadata,
            "force_checks": force_checks,
        }
        output_path = Path(args.probe_json).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(probe_record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"Probe JSON: {output_path}")
    if failed:
        print(
            f"FAILED: at least one force difference exceeds "
            f"{args.fd_tolerance:.3e} Eh/bohr",
            file=sys.stderr,
        )
        return 2
    return 0


def self_test() -> int:
    synthetic_output = """
 QED-HF ground-state energy : -7.612345678901D+01
 Dipole moment [Debye]:
 Conversion factor from Debye a.u.: 0.393430270
 Comp.         Electronic           Nuclear             Total
 x         -1.0000000000      2.0000000000      1.0000000000
 y         -2.0000000000      4.0000000000      2.0000000000
 z         -3.0000000000      6.0000000000      3.0000000000
 eT terminated successfully!
"""
    assert math.isclose(parse_energy(synthetic_output), -76.12345678901)
    assert parse_total_dipole_debye(synthetic_output) == (1.0, 2.0, 3.0)
    template = (
        "- boson\n  coherent state: {@@COHERENT_Z@@}\n"
        "  frequency: {@@FREQUENCY@@}\n  coupling: {@@COUPLING@@}\n"
        "  polarization: {@@POL_X@@, @@POL_Y@@, @@POL_Z@@}\n"
        "- geometry\n@@GEOMETRY_BOHR@@\n"
    )
    rendered = render_input(
        template,
        ["H"],
        [0.0, 0.0, 0.0],
        0.5,
        0.05,
        [0.0, 0.0, 1.0],
        0.3,
    )
    assert "coherent state: {0.15}" in rendered
    restart_template = "- solver scf\n  max iterations: 20\n- do\n  hf\n"
    assert "\n  restart\n" in enable_scf_restart(restart_template)
    with tempfile.TemporaryDirectory() as directory:
        gradient = Path(directory) / "eT.molecular_gradient.out"
        gradient.write_text(
            "atom       x                 y                 z\n"
            "O  1.0D-3  2.0D-3  3.0D-3\n"
            "H -1.0D-3 -1.0D-3 -1.5D-3\n"
            "H  0.0D+0 -1.0D-3 -1.5D-3\n",
            encoding="utf-8",
        )
        values = parse_gradient_file(gradient, 3)
        assert values == [
            1.0e-3,
            2.0e-3,
            3.0e-3,
            -1.0e-3,
            -1.0e-3,
            -1.5e-3,
            0.0,
            -1.0e-3,
            -1.5e-3,
        ]
    omega = 0.5
    coupling = 0.05
    q = 0.3
    mu = 0.9
    force = -(omega**2 * q + omega * coupling * mu)
    assert math.isclose(force, -0.0975)
    print("Self-test passed")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    connection = parser.add_mutually_exclusive_group()
    connection.add_argument(
        "--unix",
        metavar="NAME_OR_PATH",
        default="h2o_cbohf",
        help="UNIX socket name (mapped to /tmp/ipi_NAME) or absolute path",
    )
    connection.add_argument("--host", help="TCP host instead of a UNIX socket")
    parser.add_argument("--port", type=int, default=31415, help="TCP port")
    parser.add_argument("--connect-timeout", type=float, default=120.0)
    parser.add_argument("--et-launch", default=os.environ.get("ET_LAUNCH", "eT_launch.py"))
    parser.add_argument("--et-binary", help="Run the compiled eT directly (optional)")
    parser.add_argument("--basis-dir", help="eT ao_basis directory for direct execution")
    parser.add_argument("--template", default=str(Path(__file__).resolve().parent / "templates/et_cbohf_template.inp"))
    parser.add_argument("--workdir", default="et_runs")
    parser.add_argument(
        "--nuclear-symbols", nargs="+", default=["O", "H", "H"],
        help=(
            "real nuclear symbols in i-PI order; for --probe, the single value "
            "auto infers all symbols preceding the final L particle"
        ),
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--frequency", type=float, default=0.01666)
    parser.add_argument("--coupling", type=float, default=0.005)
    parser.add_argument(
        "--polarization", type=float, nargs=3, default=[0.0, 1.0, 0.0]
    )
    parser.add_argument(
        "--photon-axis",
        choices=sorted(AXES),
        default="x",
        help="Cartesian coordinate of pseudo-atom L used to store q",
    )
    parser.add_argument(
        "--gradient-mode",
        choices=["analytic", "finite-difference"],
        default="analytic",
        help="analytic requires the included eT fixed-displacement patch",
    )
    parser.add_argument(
        "--scf-restart",
        choices=["auto", "off"],
        default="auto",
        help=(
            "auto reuses the last successful per-client SCF orbitals and retries "
            "from SAD on failure; off always starts from SAD"
        ),
    )
    parser.add_argument(
        "--keep-every",
        type=int,
        default=0,
        metavar="N",
        help=(
            "retain a full calculation snapshot every N completed force "
            "evaluations; 0 keeps only the rolling current calculation"
        ),
    )
    parser.add_argument(
        "--no-failure-snapshots",
        action="store_true",
        help="do not preserve failed restarted-SCF/SAD attempts",
    )
    parser.add_argument(
        "--fd-nuclear-step",
        type=float,
        default=1.0e-4,
        help="nuclear step for --gradient-mode finite-difference (bohr)",
    )
    parser.add_argument(
        "--gradient-glob",
        default="*grad*",
        help="recursive glob selecting eT's separate gradient output",
    )
    parser.add_argument(
        "--energy-regex",
        help="override energy parser; regex must contain one capturing group",
    )
    parser.add_argument("--probe", metavar="XYZ", help="run one eT evaluation, no socket")
    parser.add_argument(
        "--probe-q",
        type=float,
        help="with --probe, override the photon coordinate q (atomic units)",
    )
    parser.add_argument(
        "--probe-json",
        help="with --probe, write energy, forces, dipoles, timing, and checks as JSON",
    )
    parser.add_argument(
        "--check-force",
        metavar="ATOM,AXIS",
        action="append",
        default=[],
        help="with --probe, compare a reported force to a central difference; repeatable",
    )
    parser.add_argument("--fd-check-step", type=float, default=1.0e-4)
    parser.add_argument("--fd-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--self-test", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.self_test:
        return self_test()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    if args.frequency <= 0.0 or args.coupling < 0.0:
        raise ValueError("Frequency must be positive and coupling non-negative")
    if args.fd_nuclear_step <= 0.0 or args.fd_check_step <= 0.0:
        raise ValueError("Finite-difference steps must be positive")
    if args.keep_every < 0:
        raise ValueError("--keep-every must be non-negative")
    if args.nuclear_symbols == ["auto"]:
        if not args.probe:
            raise ValueError("--nuclear-symbols auto is supported only with --probe")
        probe_symbols, _ = parse_xyz(Path(args.probe))
        if not probe_symbols or probe_symbols[-1] != "L" or "L" in probe_symbols[:-1]:
            raise ValueError("Auto-detected probe must contain one final L particle")
        args.nuclear_symbols = probe_symbols[:-1]
    evaluator = ETEvaluator(args, args.nuclear_symbols)
    if args.probe:
        return run_probe(args, evaluator)
    if args.probe_q is not None:
        raise ValueError("--probe-q requires --probe")
    if args.check_force:
        raise ValueError("--check-force requires --probe")
    sock = (
        connect_inet(args.host, args.port, args.connect_timeout)
        if args.host
        else connect_unix(args.unix, args.connect_timeout)
    )
    try:
        serve(sock, evaluator)
    finally:
        sock.close()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
