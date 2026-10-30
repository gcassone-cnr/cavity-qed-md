#!/usr/bin/env python3
"""Tiny eT-shaped executable used only by the interface integration test."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--omp")
    parser.add_argument("-of", required=True)
    parser.add_argument("-nt", action="store_true")
    parser.add_argument("-s", action="store_true")
    parser.add_argument("-save")
    parser.add_argument("-load")
    args = parser.parse_args()
    input_text = Path(args.input).read_text(encoding="utf-8")
    if (Path.cwd().parent / "always_fail").exists():
        print("deliberate unconditional eT failure", file=sys.stderr)
        raise SystemExit(10)
    if args.load:
        fail_once = Path(args.load) / "fail_on_load_once"
        if fail_once.exists():
            fail_once.unlink()
            print("deliberate restarted-SCF failure", file=sys.stderr)
            raise SystemExit(9)
    iterations = 2 if args.load else 13
    Path(args.of).write_text(
        f"""
 Total energy: -76.000000000000
 Convergence criterion met in {iterations} iterations!
 Dipole moment [Debye]:
 Conversion factor from Debye a.u.: 0.393430270
 Comp.         Electronic           Nuclear             Total
 x          0.0000000000      0.0000000000      0.0000000000
 y          1.0000000000      1.0000000000      2.0000000000
 z          0.0000000000      0.0000000000      0.0000000000
 eT terminated successfully!
""",
        encoding="utf-8",
    )
    if "molecular gradient" in input_text.lower():
        Path("eT.molecular_gradient.out").write_text(
            """
 - Molecular gradient (Hartree/bohr)
      Center                 x                y                z
 ----------------------------------------------------------------------
     1 O          0.010000000000   0.020000000000   0.030000000000
     2 H         -0.005000000000  -0.010000000000  -0.015000000000
     3 H         -0.005000000000  -0.010000000000  -0.015000000000
 ----------------------------------------------------------------------
""",
            encoding="utf-8",
        )
    if args.save:
        restart_directory = Path(args.save) / f"{Path(args.of).stem}_restart_files"
        restart_directory.mkdir(parents=True)
        (restart_directory / "orbital_coefficients_qed-rhf").write_text(
            "fake orbitals\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()
