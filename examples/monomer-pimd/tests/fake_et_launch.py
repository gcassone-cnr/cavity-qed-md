#!/usr/bin/env python3
"""Tiny eT-shaped executable used only by the interface integration test."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--omp")
    parser.add_argument("-of", required=True)
    args = parser.parse_args()
    Path(args.input).read_text(encoding="utf-8")
    Path(args.of).write_text(
        """
 Total energy: -76.000000000000
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


if __name__ == "__main__":
    main()
