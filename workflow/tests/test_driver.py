#!/usr/bin/env python3
"""Unit tests that do not require eT or i-PI installations."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import socket
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest


PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "et_cbohf_ipi_driver", PROJECT / "et_cbohf_ipi_driver.py"
)
assert SPEC is not None and SPEC.loader is not None
driver = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = driver
SPEC.loader.exec_module(driver)


def arguments(workdir: str, gradient_mode: str = "analytic") -> argparse.Namespace:
    return argparse.Namespace(
        et_launch="unused",
        template=str(PROJECT / "templates" / "et_cbohf_template.inp"),
        workdir=workdir,
        frequency=0.5,
        coupling=0.05,
        polarization=[0.0, 1.0, 0.0],
        photon_axis="x",
        threads=1,
        energy_regex=None,
        gradient_glob="*grad*",
        gradient_mode=gradient_mode,
        fd_nuclear_step=1.0e-4,
        scf_restart="auto",
        keep_every=0,
        no_failure_snapshots=False,
    )


class DriverTests(unittest.TestCase):
    def test_probe_q_option(self) -> None:
        args = driver.build_parser().parse_args(
            ["--probe", "geometry.xyz", "--probe-q", "0.5"]
        )
        self.assertEqual(args.probe_q, 0.5)

    def test_probe_json_contains_forces_and_timing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            xyz = Path(directory) / "probe.xyz"
            xyz.write_text(
                "4\nprobe\n"
                "O 0.0 0.0 0.0\n"
                "H 0.0 0.0 1.0\n"
                "H 0.0 1.0 0.0\n"
                "L 0.0 0.0 0.0\n",
                encoding="utf-8",
            )
            output = Path(directory) / "probe.json"
            args = driver.build_parser().parse_args(
                ["--probe", str(xyz), "--probe-json", str(output)]
            )
            evaluator = SimpleNamespace(
                nuclear_symbols=["O", "H", "H"],
                photon_axis=0,
                evaluate=lambda positions: (
                    -76.0,
                    [0.0] * len(positions),
                    {
                        "nuclear_gradient_mode": "analytic",
                        "q_au": 0.0,
                        "polarized_dipole_au": 1.0,
                        "et_launcher_wall_seconds": 2.5,
                    },
                ),
            )
            self.assertEqual(driver.run_probe(args, evaluator), 0)
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["symbols"], ["O", "H", "H", "L"])
            self.assertEqual(len(payload["forces_au"]), 12)
            self.assertEqual(payload["metadata"]["et_launcher_wall_seconds"], 2.5)

    def test_render_converts_q_to_coherent_z(self) -> None:
        text = driver.render_input(
            (PROJECT / "templates" / "et_cbohf_template.inp").read_text(encoding="utf-8"),
            ["O", "H", "H"],
            [0.0] * 9,
            0.5,
            0.05,
            [0.0, 1.0, 0.0],
            0.3,
        )
        self.assertIn("coherent state: {0.15}", text)
        self.assertNotIn("@@", text)

    def test_render_enables_restart_and_can_omit_gradient(self) -> None:
        text = driver.render_input(
            (PROJECT / "templates" / "et_cbohf_template.inp").read_text(encoding="utf-8"),
            ["O", "H", "H"],
            [0.0] * 9,
            0.5,
            0.05,
            [0.0, 1.0, 0.0],
            0.3,
            include_gradient=False,
            scf_restart=True,
        )
        self.assertIn("- solver scf\n  restart\n", text)
        self.assertNotIn("molecular gradient", text.lower())

    def test_dipole_parser_uses_total_column(self) -> None:
        text = """
 Dipole moment [Debye]:
 Comp.         Electronic           Nuclear             Total
 x         -4.0       5.0       1.0
 y         -5.0       7.0       2.0
 z         -6.0       9.0       3.0
 """
        self.assertEqual(driver.parse_total_dipole_debye(text), (1.0, 2.0, 3.0))

    def test_analytic_joint_force(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            evaluator = driver.ETEvaluator(arguments(directory), ["O", "H", "H"])

            def fake_run_et(
                calculation_directory: Path,
                nuclear_positions_bohr: list[float],
                q: float,
                need_gradient: bool,
            ) -> object:
                del calculation_directory
                energy = (
                    0.25 * sum(value * value for value in nuclear_positions_bohr)
                    + 0.5 * evaluator.frequency**2 * q**2
                    + evaluator.frequency * evaluator.coupling * 2.0 * q
                )
                gradient = (
                    [0.5 * value for value in nuclear_positions_bohr]
                    if need_gradient
                    else None
                )
                return driver.ETPoint(
                    energy=energy,
                    dipole_au=(0.0, 2.0, 0.0),
                    gradient=gradient,
                    gradient_file="gradient.out" if need_gradient else None,
                )

            evaluator._run_et = fake_run_et
            positions = [0.1] * 9 + [0.3, 0.0, 0.0]
            energy, forces, metadata = evaluator.evaluate(positions)
            self.assertTrue(math.isfinite(energy))
            self.assertEqual(len(forces), 12)
            self.assertAlmostEqual(forces[0], -0.05)
            self.assertAlmostEqual(forces[9], -0.125)
            self.assertEqual(forces[10:], [0.0, 0.0])
            self.assertAlmostEqual(metadata["polarized_dipole_au"], 2.0)

    def test_finite_difference_nuclear_gradient(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            evaluator = driver.ETEvaluator(
                arguments(directory, "finite-difference"), ["O", "H", "H"]
            )

            def fake_run_et(
                calculation_directory: Path,
                nuclear_positions_bohr: list[float],
                q: float,
                need_gradient: bool,
            ) -> object:
                del calculation_directory, need_gradient
                energy = 0.25 * sum(
                    value * value for value in nuclear_positions_bohr
                )
                return driver.ETPoint(
                    energy=energy,
                    dipole_au=(0.0, 0.0, 0.0),
                    gradient=None,
                    gradient_file=None,
                )

            evaluator._run_et = fake_run_et
            positions = [0.2] * 9 + [0.0, 0.0, 0.0]
            _, forces, _ = evaluator.evaluate(positions)
            for force in forces[:9]:
                self.assertAlmostEqual(force, -0.1, places=10)

    def test_force_payload(self) -> None:
        left, right = socket.socketpair()
        try:
            driver.send_force(left, -1.5, [0.1, 0.2, 0.3], {"q": 0.4})
            self.assertEqual(driver.recv_exact(right, 12), driver.header("FORCEREADY"))
            self.assertAlmostEqual(
                struct.unpack("=d", driver.recv_exact(right, 8))[0], -1.5
            )
            self.assertEqual(
                struct.unpack("=i", driver.recv_exact(right, 4))[0], 1
            )
            self.assertEqual(
                struct.unpack("=3d", driver.recv_exact(right, 24)),
                (0.1, 0.2, 0.3),
            )
            driver.recv_exact(right, 72)
            extra_size = struct.unpack("=i", driver.recv_exact(right, 4))[0]
            extra = json.loads(driver.recv_exact(right, extra_size))
            self.assertEqual(extra, {"q": 0.4})
        finally:
            left.close()
            right.close()

    def test_fake_et_subprocess_and_parsers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = arguments(directory)
            args.et_launch = str(PROJECT / "tests" / "fake_et_launch.py")
            evaluator = driver.ETEvaluator(args, ["O", "H", "H"])
            energy, forces, metadata = evaluator.evaluate(
                [0.0] * 9 + [0.0, 0.0, 0.0]
            )
            self.assertEqual(energy, -76.0)
            self.assertAlmostEqual(forces[0], -0.01)
            self.assertAlmostEqual(forces[1], -0.02)
            self.assertAlmostEqual(
                metadata["polarized_dipole_au"],
                2.0 * driver.DEBYE_TO_AU,
            )
            self.assertAlmostEqual(
                forces[9],
                -0.5 * 0.05 * 2.0 * driver.DEBYE_TO_AU,
            )
            self.assertFalse(metadata["scf_restart_loaded"])
            self.assertEqual(metadata["scf_iterations"], 13)
            _, _, second = evaluator.evaluate([0.0] * 9 + [0.0, 0.0, 0.0])
            self.assertTrue(second["scf_restart_loaded"])
            self.assertFalse(second["scf_restart_fallback"])
            self.assertEqual(second["scf_iterations"], 2)
            self.assertTrue((Path(directory) / "current" / "et.out").is_file())
            self.assertTrue(
                (
                    Path(directory)
                    / "scf_restart"
                    / "orbital_coefficients_qed-rhf"
                ).is_file()
            )
            self.assertEqual(list(Path(directory).glob("step_*")), [])

    def test_restart_failure_falls_back_to_sad_and_keeps_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = arguments(directory)
            args.et_launch = str(PROJECT / "tests" / "fake_et_launch.py")
            evaluator = driver.ETEvaluator(args, ["O", "H", "H"])
            positions = [0.0] * 9 + [0.0, 0.0, 0.0]
            evaluator.evaluate(positions)
            (Path(directory) / "scf_restart" / "fail_on_load_once").touch()
            _, _, metadata = evaluator.evaluate(positions)
            self.assertFalse(metadata["scf_restart_loaded"])
            self.assertTrue(metadata["scf_restart_fallback"])
            self.assertEqual(metadata["scf_iterations"], 13)
            failures = list((Path(directory) / "failures").glob("*restart_failed"))
            self.assertEqual(len(failures), 1)
            self.assertTrue((failures[0] / "calculation" / "launcher.stderr").is_file())

    def test_double_failure_preserves_last_known_good_restart(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = arguments(directory)
            args.et_launch = str(PROJECT / "tests" / "fake_et_launch.py")
            evaluator = driver.ETEvaluator(args, ["O", "H", "H"])
            positions = [0.0] * 9 + [0.0, 0.0, 0.0]
            evaluator.evaluate(positions)
            orbitals = (
                Path(directory)
                / "scf_restart"
                / "orbital_coefficients_qed-rhf"
            )
            before = orbitals.read_bytes()
            (Path(directory) / "always_fail").touch()
            with self.assertRaises(RuntimeError):
                evaluator.evaluate(positions)
            self.assertEqual(orbitals.read_bytes(), before)
            failures = list((Path(directory) / "failures").iterdir())
            self.assertEqual(len(failures), 2)

    def test_sparse_snapshot_interval(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = arguments(directory)
            args.keep_every = 2
            args.et_launch = str(PROJECT / "tests" / "fake_et_launch.py")
            evaluator = driver.ETEvaluator(args, ["O", "H", "H"])
            positions = [0.0] * 9 + [0.0, 0.0, 0.0]
            evaluator.evaluate(positions)
            self.assertFalse((Path(directory) / "snapshots").exists())
            evaluator.evaluate(positions)
            self.assertTrue(
                (Path(directory) / "snapshots" / "force_000002" / "et.out").is_file()
            )

    def test_energy_only_evaluation_omits_gradient_task(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            args = arguments(directory)
            args.et_launch = str(PROJECT / "tests" / "fake_et_launch.py")
            evaluator = driver.ETEvaluator(args, ["O", "H", "H"])
            energy = evaluator.evaluate_energy([0.0] * 9 + [0.0, 0.0, 0.0])
            self.assertEqual(energy, -76.0)
            current = Path(directory) / "current"
            self.assertNotIn(
                "molecular gradient",
                (current / "et.inp").read_text(encoding="utf-8").lower(),
            )
            self.assertFalse((current / "eT.molecular_gradient.out").exists())


if __name__ == "__main__":
    unittest.main()
