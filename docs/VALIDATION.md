# Validation

## Electronic force consistency

Compare the analytic nuclear and photon forces with central differences of the relaxed CBO-HF energy. Run the checks with the compiled eT executable and the electronic template used for dynamics.

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

The independent molecular calculation for the v4 fixed-displacement gradient found maximum nuclear force differences below 1e-9 a.u. at finite coupling for q = 0 and 0.5. See the [gradient correction and reference calculation](GRADIENT_FIX_V4.md) for the comparison and its scope.

## Interface and build checks

| Check | Coverage |
|---|---|
| Driver tests | Input staging, energy and force parsing, socket protocol, SCF continuation and cold-start recovery |
| i-PI integration | Monomer and cluster dynamics with harmonic test forces, including the one-bead Eckart constraints |
| eT source generation | Generation from the bundled eT archive, repeat application and preservation of the patched QED-HF source |
| Distribution integrity | SHA-256 checksums, archive contents and local website links |

The GitHub Actions workflow runs these checks. The harmonic integration tests exercise the interface and integrator. Use the force suite above to validate the compiled electronic engine. Historical workstation records are in `validation/`. [GitHub Actions](https://github.com/gcassone-cnr/cavity-qed-md/actions)

## Simulation convergence

The examples cover one explicit cavity mode and finite molecular systems on a CBO-HF surface. Compare NVE energy conservation at matched durations with at least two timesteps. For equilibrium sampling, converge the basis, SCF tolerance, bead count and sampling length, and report statistical uncertainties. The monomer PILE-L input provides a short joint nuclear–photon PIMD run. The accompanying [arXiv article](https://arxiv.org/abs/2611.09667) describes the scientific method and benchmarks.
