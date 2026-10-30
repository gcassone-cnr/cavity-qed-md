#!/usr/bin/env python3
"""Launch i-PI with compatibility for its Eckart constrained-MD API.

i-PI releases 3.1.12, 3.2.0 and 3.3.0 contain stale interfaces in the
constrained-dynamics path: (1) ``fixatoms`` versus ``fixatoms_dof`` and (2) the
removed ``forces_mts(level)`` method versus ``mts_forces[level].f``. Ordinary
dynamics is unaffected. This wrapper applies the minimal in-memory
compatibility aliases and then executes the installed ``i-pi`` command without
modifying the installation.
"""

from __future__ import annotations

import inspect
import os
import runpy
import shutil
import sys


def apply_compatibility() -> list[str]:
    from ipi.engine.forces import Forces
    from ipi.engine.motion.constrained_dynamics import (
        ConstrainedDynamics,
        NVEConstrainedIntegrator,
    )
    from ipi.engine.motion.motion import Motion
    from ipi.utils.depend import dstrip

    applied = []
    cd_parameters = inspect.signature(ConstrainedDynamics.__init__).parameters
    motion_parameters = inspect.signature(Motion.__init__).parameters
    cd_needs_alias = "fixatoms_dof" not in cd_parameters and "fixatoms" in cd_parameters
    motion_needs_alias = (
        "fixatoms" not in motion_parameters and "fixatoms_dof" in motion_parameters
    )

    if cd_needs_alias:
        original_cd_init = ConstrainedDynamics.__init__

        if motion_needs_alias:
            original_motion_init = Motion.__init__

            def motion_init_compat(
                self, fixcom=False, fixatoms_dof=None, fixatoms=None
            ):
                if fixatoms_dof is None:
                    fixatoms_dof = fixatoms
                return original_motion_init(
                    self, fixcom=fixcom, fixatoms_dof=fixatoms_dof
                )

            Motion.__init__ = motion_init_compat

        def constrained_init_compat(
            self, *args, fixatoms_dof=None, fixatoms=None, **kwargs
        ):
            if fixatoms is None:
                fixatoms = fixatoms_dof
            return original_cd_init(self, *args, fixatoms=fixatoms, **kwargs)

        ConstrainedDynamics.__init__ = constrained_init_compat
        applied.append("fixed-DOF parameter alias")

    if not hasattr(Forces, "forces_mts"):
        def constrained_step_b_compat(self, level=0):
            force = dstrip(self.forces.mts_forces[level].f)[0]
            self.beads.p[0] += force * self.pdt[level]

        NVEConstrainedIntegrator.step_B = constrained_step_b_compat
        applied.append("MTS force accessor")

    return applied


def find_real_ipi() -> str:
    wrapper = os.path.realpath(sys.argv[0])
    candidates = []
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        if not directory:
            continue
        candidate = os.path.realpath(os.path.join(directory, "i-pi"))
        if candidate != wrapper and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            candidates.append(candidate)
    if not candidates:
        fallback = shutil.which("i-pi")
        if fallback and os.path.realpath(fallback) != wrapper:
            candidates.append(os.path.realpath(fallback))
    if not candidates:
        raise RuntimeError(
            "Cannot find the real i-pi executable. Activate the i-PI virtual "
            "environment before running this wrapper."
        )
    return candidates[0]


def main() -> None:
    patches = apply_compatibility()
    if patches:
        print(
            "i-pi-eckart-compat: applied constrained-dynamics compatibility: "
            + ", ".join(patches),
            file=sys.stderr,
        )
    ipi_executable = find_real_ipi()
    sys.argv[0] = ipi_executable
    runpy.run_path(ipi_executable, run_name="__main__")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"i-pi-eckart-compat: ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
