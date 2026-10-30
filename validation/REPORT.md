**Validation report — 18 September 2026**

The regression suite ran 21 tests: **20 passed, 1 skipped**. The skipped test is the complete native-driver/i-PI UNIX-socket integration test; this execution environment denies UNIX socket creation. It is included in `workflow/tests/test_workstation.py` and will run on a normal Linux workstation during installation. A skipped test is not a pass.

Tests ran with Python 3.12.14, i-PI 3.2.0, NumPy 2.2.6 and SciPy 1.15.3. Full interpreter details are in `environment.json`; the test transcript is `tests.log`.

**What was exercised**

- Original restart, failure-preservation, snapshot-retention and probe JSON regression cases.
- Direct native-binary input/output handling with a harmonic electronic test double.
- Nuclear and photon analytic forces against central energy differences for that test double.
- Cold starts, warm starts and a deliberately corrupted restart followed by cold fallback.
- Rejection of unsupported DFT templates in the CBO driver.
- Normal idle socket EOF handling and rejection of truncated protocol messages.
- Correct thermal energy units, mass-weighted momentum projection and retained photon momentum.
- Preservation of the original Eckart reference for a deformed checkpoint configuration.
- Rejection of incompatible fixed Cartesian DOFs in Eckart input.
- Real pinned i-PI constrained integration with its internal harmonic force field, for the monomer and the 21-water cluster.

The two real i-PI harmonic integrations each used 30 steps at 0.1 fs. These are integrator checks with a harmonic potential, **not electronic CBO-HF calculations**:

| System | Final Eckart position-residual norm | Final tangent-velocity residual norm | Conserved-energy peak-to-peak, eV |
|---|---:|---:|---:|
| One water + L | 3.35e-16 | 1.08e-20 | 3.05e-5 |
| 21 waters + L | 1.31e-14 | 7.25e-20 | 2.10e-6 |

Residuals use the normalization and internal atomic units of i-PI's constraint implementation. Unused photon y/z positions and momenta remained exactly zero. Energy ranges are calculated from the printed properties output and inherit its precision. Raw metric JSON files and i-PI logs are included.

The driver self-test passed. Python compilation and Bash syntax checks passed. The recovered gradient patch applies with zero fuzz to the bundled eT source in a dry run. A launcher dry run identified 63 nuclei, the photon and global Eckart constraints correctly (`staging.log`).

**Not established here**

- Compilation/linking of eT, libcint and LibXC on the target workstation. A Fortran compiler was unavailable locally, and system-package installation was restricted.
- Numerical electronic-force accuracy of the recovered eT gradient patch in the selected basis, geometry, q and coupling regime.
- Agreement of the real eT direct and launcher execution paths. `compare_backends.py` provides the local gate.
- Full UNIX socket integration, ordinary-NVT end-to-end execution, sustained production stability, IR peak convergence or actual 21-water throughput.
- CBO-DFT dynamics or forces. The included PBE input is an ordinary DFT single-point example.

The installer and README provide the remaining workstation checks. A patch installer saying “already applied” only identifies source state; it does not certify the compiled executable or its gradients. Mandatory sparse force checks are supplemented by `validate_forces.py --all-nuclear`, selected cluster snapshot checks and NVE timestep comparisons. Numerical speedup is intentionally left to measurements on the new machine.
