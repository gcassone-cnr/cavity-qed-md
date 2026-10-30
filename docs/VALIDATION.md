# Validation status — 30 October 2026

On 30 October 2026, Giuseppe Cassone confirmed the corrected v4 installation and documented test workflow as working on his workstation. This is a maintainer report; raw logs from that successful target-machine run are not bundled. The separate packaging and independent molecular checks below retain their own scope. See `validation/WORKSTATION_STATUS.md`.

## Checks rerun during preparation

| Suite | Result | What it establishes |
|---|---|---|
| Workstation suite (`workflow/tests`) | 21 tests discovered: 20 passed, 1 skipped | Parsers, energy/force agreement for a harmonic electronic test double, SCF restart/fallback, input staging, and real i-PI internal-harmonic Eckart integration |
| Original monomer/PIMD driver suite | 7 tests passed | Protocol payloads, parsing, q-to-coherent mapping and analytic/finite-difference agreement for the test model |
| Full eT/i-PI UNIX-socket test | Skipped by the workstation suite | This execution environment prohibits UNIX socket creation; the test is included for the target workstation |

The real i-PI internal-harmonic test is an integration test of i-PI, not a real eT electronic-force calculation. The regenerated report does not turn earlier skipped tests into passes.

## Distribution and website checks

The original driver, templates, examples and vendored eT bytes were compared with their recovered archives and found unchanged. The gradient patch and its application script have since been corrected to v4 as described below. Source checksums, the downloadable archive, local website resources, Python syntax and Bash syntax passed. The landing page was visually inspected in a desktop browser; its MD/PIMD/analysis tabs and copy-button feedback were checked. Mobile styling is included but was not separately inspected at a mobile viewport.

The reported installer failure at missing `src/eri_tool/generated_eri_tool/eri_tool_c_class.F90` was traced to the omitted eT source-generation step. The corrected installer calls eT's own `autogenerate(root)` before CMake. `scripts/check_et_generation.py` extracts a clean vendored snapshot, applies the QED-HF patch, executes the installer's generation block twice, and checks every explicit Fortran source in the main and libcint CMake source lists. It also verifies that regeneration preserves the patched QED-HF source. These generation checks do not establish successful compilation or real electronic forces.

## Required checks with the installed eT executable

The first reported real-eT test passed lambda = 0 but failed nuclear forces at lambda = 0.01, q = 0, while the photon force agreed to about 2.2e-9 a.u. An independent molecular CBO-HF calculation reproduced the reported energy and old force-error pattern. The corrected fixed-displacement DSE gradient gives maximum nuclear errors below 1e-9 a.u. at q = 0 and 0.5 in that independent check. Zero-coupling checks also passed. See [v4 correction](GRADIENT_FIX_V4.md) and `validation/independent_cbohf_v4.log`. The maintainer subsequently reported the corrected v4 installation and test workflow working. Every new installation must repeat the real-eT acceptance checks; harmonic test-double passes alone do not validate the Fortran correction.

From the installed `workflow/` directory:

```bash
python3 -m unittest discover -s tests -v
python3 compare_backends.py
python3 validate_forces.py --all-nuclear

python3 validate_forces.py --xyz examples/water21.xyz \
  --threads 16 --frequency 0.0187505 --couplings 0.01 \
  --output force_validation_21
```

The monomer force suite includes zero coupling and nonzero q. Its default absolute force tolerance is 1e-5 atomic units. The cluster command checks selected components unless `--all-nuclear` is added. Expand configurations, q values and coupling ranges to cover the actual application. Use the matching executable, electronic template and cavity parameters.

Archive the compiler/build identity, `installation.json`, the force-validation outputs and representative NVE energy-conservation results with a release. Compare at least two timesteps at matched physical duration; converge bead count, basis, SCF settings, sampling length and uncertainties for scientific claims.

## Scope of the distributed examples

- Ordinary monomer/21-water MD inputs and one-bead global Eckart NVE are included.
- The original monomer joint nuclear–photon PIMD input is included; a fresh end-to-end PIMD/eT run was not performed here.
- The user reports a PLUMED interface, but its current source and runnable example were not recovered for this release.
- CBO-DFT/correlated-method dynamics, converged reaction calculations, and vacuum-induced catalytic effects are not established by this package.

The original 18 September validation report and logs remain in `validation/`. Read their date and test-model descriptions when comparing them with this report.
