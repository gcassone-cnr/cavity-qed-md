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
        template=str(PROJECT / "et_cbohf_template.inp"),
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
    )


class DriverTests(unittest.TestCase):
    def test_probe_q_option(self) -> None:
        args = driver.build_parser().parse_args(
            ["--probe", "geometry.xyz", "--probe-q", "0.5"]
        )
        self.assertEqual(args.probe_q, 0.5)

    def test_render_converts_q_to_coherent_z(self) -> None:
        text = driver.render_input(
            (PROJECT / "et_cbohf_template.inp").read_text(encoding="utf-8"),
            ["O", "H", "H"],
            [0.0] * 9,
            0.5,
            0.05,
            [0.0, 1.0, 0.0],
            0.3,
        )
        self.assertIn("coherent state: {0.15}", text)
        self.assertNotIn("@@", text)

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


if __name__ == "__main__":
    unittest.main()
