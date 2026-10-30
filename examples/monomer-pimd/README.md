# Joint nuclear–photon CBO-HF dynamics for H2O

This project couples a fixed-displacement eT QED-HF calculation to an ordinary
i-PI socket. It propagates the three water nuclei and one explicit cavity
coordinate \(q\). The same interface can run classical NVE dynamics or joint
nuclear–photon PIMD.

## What the source audit found

eT already reads a fixed cavity coordinate (`boson / cbo`) and converts it to a
coherent displacement,

\[
z=q\sqrt{\omega/2}.
\]

However, in the uploaded development source the ordinary CBO-HF gradient path
does not differentiate all cavity terms. The QED-HF gradient contains the
dipole-self-energy derivative, but lacks the explicit derivatives of the
fixed coherent-state terms. Using either stock analytic gradient for
fixed-\(q\) MD would therefore be unsafe.

This project uses `method qed-hf` plus a fixed `coherent state`, which avoids
eT's separate CBO embedding path and retains the complete one- and two-electron
dipole self-energy terms. The included small patch adds the missing
fixed-displacement terms to the QED-HF molecular gradient. The fixed-displacement branch includes both the coherent interaction and the
DSE Hartree derivative, contracted with eT's complete polarized-dipole derivative
tensor including the direct nuclear-dipole derivative. The optimized-displacement
overlap term is kept only for optimized bosons; fixed displacement uses the
actual harmonic overlap term. These v4 corrections are needed even at q=0
when the coupling is nonzero. The optimized-boson branch is unchanged.

The patch is source-reviewed and the Python interface is unit-tested here, but
this environment has no Fortran compiler or eT dependencies. The patched
Fortran has therefore **not** been compiled or numerically validated here.
`run.sh` makes both a nuclear and photon finite-difference force check mandatory
before it launches i-PI.

## Model and coordinate convention

For one cavity mode, eT's sign convention gives

\[
E(R,q)=E_\mathrm{el}(R;q)
       +\frac12\omega^2q^2
       +\omega\lambda q\,\langle\boldsymbol{\mu}\cdot\mathbf e\rangle
       +E_\mathrm{DSE},
\]

and therefore

\[
F_q=-\omega^2q-\omega\lambda
       \langle\boldsymbol{\mu}\cdot\mathbf e\rangle .
\]

The Python client obtains the total molecular dipole from eT and evaluates this
photon force analytically. It sends only O/H/H to eT. The fourth i-PI particle,
labelled `L`, represents the photon:

- `L.x` is \(q\);
- `L.y` and `L.z` are frozen;
- the physical cavity polarization remains the y direction, along water's
  permanent-dipole/C2 axis;
- the ordinary `ffsocket` is used. Do not add `ffcavphsocket`, because that
  would add a second cavity model.

The supplied mode has `omega = 0.01666 Eh` (about 3656 cm-1),
`lambda = 0.005 au`, and QED-HF/aug-cc-pVDZ.

## 1. Patch and rebuild eT

Keep an unmodified copy or a git branch, then run:

```bash
./apply_et_patch.sh /absolute/path/to/eT-development
cd /absolute/path/to/eT-development
cmake --build build-gcc -j 8
```

If the source was already configured, only the final build command should be
needed. The patch touches:

```text
src/wavefunctions/qed_hf/qed_hf_class.F90
```

It is also available directly as
`patches/et-fixed-coherent-gradient.patch`.

If an earlier version of this package was already patched into eT,
`apply_et_patch.sh` detects and upgrades it automatically.

## 2. Run the mandatory force check and 10 fs NVE test

```bash
export ET_LAUNCH=/absolute/path/to/eT-development/build-gcc/eT_launch.py
export ET_THREADS=4
./run.sh
```

Each launch creates `runs/<timestamp>/`. Before i-PI starts, `run.sh` checks
the oxygen y force and photon `L.x` force against central energy differences
at both \(q=0\) and \(q=0.5\) au. The nonzero-\(q\) check is essential: it
validates the additional fixed-coherent-state nuclear-gradient terms.

The default tolerance is `1e-4 Eh/bohr`. If either check fails, MD does not
start. Inspect `force_check_q0.log`, `force_check_q05.log`, and the
corresponding eT calculations under `et_force_check_q0/` and
`et_force_check_q05/`.

The principal NVE outputs are:

- `h2o_cbohf.properties`;
- `h2o_cbohf.positions_au_0.xyz`;
- `ipi.log` and `et_driver_0.log`;
- `et_runs/client_0/step_XXXXXX/`, containing every eT input and output.

Trajectory coordinates are written in atomic units. In each four-particle XYZ
frame, the x coordinate of the final `L` row is \(q\).

## 3. Run eight-bead joint nuclear–photon PIMD

```bash
IPI_INPUT=ipi_cbo_pimd.xml ET_CLIENTS=8 ./run.sh
```

All nuclear coordinates and \(q\) then have eight ring-polymer beads. `L` has
the electron mass in i-PI, so its unfrozen x coordinate has the intended
canonical photon kinetic energy. `ET_CLIENTS=8` evaluates beads concurrently;
reduce it if `8 * ET_THREADS` would oversubscribe the machine.

The PIMD files are:

- `h2o_cbohf_pimd.properties`;
- `h2o_cbohf_pimd.centroid_au_0.xyz`;
- `h2o_cbohf_pimd.beads_au_0.xyz` through `_7.xyz`.

This PILE-G input samples an equilibrium quantum distribution. Its thermostat
trajectory is not a physical real-time spectrum; RPMD/TRPMD needs a separate
dynamics protocol after equilibration.

## Slow, patch-free validation mode

The driver can obtain every nuclear force by central energy differences:

```bash
ET_GRADIENT_MODE=finite-difference ./run.sh
```

For water this requires 18 additional eT calculations per force evaluation,
per bead. It is useful as a reference, but it is not practical production MD.
The fixed-displacement QED-HF energy and dipole do not require the patch.

## Controls

| Variable | Default | Meaning |
|---|---:|---|
| `ET_LAUNCH` | required | patched eT launcher |
| `IPI_INPUT` | `ipi_cbo_nve.xml` | NVE or PIMD XML |
| `ET_THREADS` | `1` | OpenMP threads per eT process |
| `ET_CLIENTS` | `1` | concurrent socket clients |
| `ET_FREQUENCY` | `0.01666` | \(\omega\), hartree |
| `ET_COUPLING` | `0.005` | \(\lambda\), atomic units |
| `ET_GRADIENT_MODE` | `analytic` | `analytic` or `finite-difference` |
| `IPI_COMMAND` | `i-pi` | i-PI executable |
| `RUN_TAG` | timestamp | output directory name |

The client has additional parser, polarization, finite-difference, TCP, and
socket options:

```bash
python3 et_cbohf_ipi_driver.py --help
```

## Validation required before scientific production

After the supplied component checks pass:

1. Check all nine nuclear force components and several nonzero \(q\) values.
2. Confirm the energy minimum in \(q\) satisfies
   \(q_\mathrm{min}=-\lambda\langle\mu_e\rangle/\omega\).
3. At `lambda=0`, recover ordinary HF nuclear forces plus an independent
   harmonic photon.
4. Check NVE conserved-energy drift at 0.1 and 0.05 fs.
5. Converge basis, SCF thresholds, PIMD bead count, trajectory length, and
   statistical uncertainty.
6. For PIMD, verify that the frozen DOFs reported by i-PI are `[10, 11]` and
   that `L.x` fluctuates while `L.y` and `L.z` remain zero.

Do not interpret a successful 10 fs smoke test as a converged cavity-QED
simulation.
