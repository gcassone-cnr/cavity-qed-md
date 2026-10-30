# Fixed-displacement nuclear-gradient correction (v4)

The original fixed-coherent-state patch was incomplete. It added the explicit
coherent-displacement derivatives but retained a QED-HF nuclear-gradient
expression specialized to an optimized photon displacement. At fixed q this
omits the dipole Hartree derivative and retains an overlap term associated with
the optimized displacement. The error is present even at q = 0 and nonzero
coupling. It is not corrected by increasing force-check tolerances.

## Evidence

For the distributed water geometry, cc-pVDZ, omega = 0.0187505 a.u. and
lambda = 0.01 a.u., an independent PySCF 2.14.0 implementation of the same
finite-basis CBO-HF energy gives -75.98957604860672 Eh at q = 0, compared with
the reported eT value -75.98957604860701 Eh. It reproduces the old nuclear
forces and their errors, including the oxygen z error of about 6.883e-5 a.u.
The photon force already agrees with the total-energy derivative.

| q (a.u.) | Maximum old nuclear force error (a.u.) | Maximum corrected error in independent check (a.u.) |
|---|---|---|
| 0 | 6.8835e-5 | 7.5e-10 |
| 0.5 | 6.8758e-5 | 8.3e-10 |

Zero-coupling nuclear and photon checks also passed, with numerical errors
below 1.3e-9 a.u. See `validation/independent_cbohf_v4.log` for the full output.
These are molecular energy/gradient checks using PySCF, not a claim that the
corrected eT Fortran executable was compiled or tested in this environment.
The cavity overlap/dipole integral derivatives in the independent check use
central differences; relaxed total-energy differences provide the reference.

## What changes in the Fortran

Let D be the idempotent RHF density (the spin-summed density is 2D),
d = lambda times the total polarized dipole AO matrix, and d_exp = 2 Tr(D d).
The fixed-displacement DSE energy contains

```math
E_{DSE}=\operatorname{Tr}(D d S^{-1}d)
       -\operatorname{Tr}(D d D d)
       +\tfrac12 d_{exp}^{2}.
```

Optimizing the coherent displacement cancels the last term against the
bilinear and harmonic contributions. Holding the displacement fixed requires
its derivative. The v4 patch restores this through the full interaction
derivative, including the direct nuclear-dipole derivative. In the code's
g = sqrt(omega/2) d convention, the coefficient of contraction with the
spin-summed density becomes `2*z + sqrt(2/omega)*d_exp`.

The `d_exp**2 / N_e` overlap contribution is retained only when the bosons are
optimized; fixed displacement uses the actual `omega*z**2` overlap term
already supplied by the earlier patch. The optimized-boson branch retains its
previous behavior. No electronic energy, photon force, basis, SCF threshold,
or validation tolerance is changed by this correction.

## Apply to an existing v3 installation

Download `et-fixed-coherent-gradient-v3-to-v4.patch` and, from the directory
containing it, run:

```bash
patch --fuzz=0 --forward -p1 \
  -d "$HOME/opt/et-ipi-release-20261030/src/eT-development" \
  < et-fixed-coherent-gradient-v3-to-v4.patch

source "$HOME/opt/et-ipi-release-20261030/env.sh"
cmake --build "$HOME/opt/et-ipi-release-20261030/build-et" --parallel 8

cd "$HOME/opt/et-ipi-release-20261030/workflow"
python3 validate_forces.py --all-nuclear --output force_validation_v4
```

The incremental patch expects v3 and refuses mismatched source hunks. For new
installations the full corrected patch is applied automatically. The updated
`workflow/apply_et_patch.sh` also recognizes existing v1/v2/v3 source and upgrades
it to v4, and recognizes an already-applied v4 patch.

The fresh validation output directory preserves the diagnostic results from
the failed v3 run. Rebuild and rerun the real-eT force checks before dynamics.

## Reproduce the independent check

In a separate Python environment with NumPy, SciPy, h5py and `pyscf==2.14.0`:

```bash
python3 scripts/check_fixed_displacement_gradient.py
```

PySCF is an optional validation dependency, not required to install or run the
eT/i-PI workflow. The generation/patch regression check also covers clean
application, upgrading from v3, and repeated v4 application.
