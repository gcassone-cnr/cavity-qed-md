# Cavity-QED Molecular Dynamics

**An eT–i-PI interface for cavity Born–Oppenheimer molecular dynamics and joint nuclear–photon path-integral sampling.**

**Start with [START_HERE.md](START_HERE.md)** for the complete consolidated v4 installation. This archive includes both the automatic eT source-generation fix and the complete gradient correction; separate repair downloads are unnecessary.

Maintained by Giuseppe Cassone. This repository distributes the corrected v4 eT–i-PI CBO-HF implementation with its matching source snapshot, examples and reproducibility checks.

**Gradient correction v4:** a reported nonzero-coupling nuclear-force mismatch was traced to missing fixed-displacement DSE derivatives. The corrected patch is included for new installations and existing v1/v2/v3 source. See [the diagnosis, independent molecular checks and upgrade commands](docs/GRADIENT_FIX_V4.md). Rebuild and validate the corrected eT executable before dynamics.

## What is included

- The fixed-coherent-state eT nuclear-gradient patch and its matching eT source snapshot.
- An i-PI socket client with analytic nuclear and photon forces, SCF restart, cold fallback, and rolling scratch.
- A workstation installer, monomer and 21-water examples, ordinary NVT/NVE preparation, and one-bead global Eckart NVE.
- The separate recovered monomer example for joint nuclear–photon PIMD.
- Force checks, backend comparisons, performance summaries, IR and photon-coordinate analysis tools.
- A self-contained [GitHub Pages website](docs/index.html) with a source download.

**Release status:** on 30 October 2026, the maintainer confirmed that the corrected v4 installation and documented test workflow work on his workstation. Interface tests and independent molecular gradient checks are included. Repeat the compiled-engine checks on each target installation; production scientific observables require their own convergence tests. See [validation](docs/VALIDATION.md).

## Install and validate the corrected v4 workflow

Target: Ubuntu 22.04/24.04 on an x86-64 CPU, Python 3.10–3.12, GNU compilers and LP64 OpenBLAS. Internet access is required for dependencies. All project sources and corrections are in this archive; manual repair patches are unnecessary.

### 1. Install system prerequisites

```bash
sudo apt update
sudo apt install -y build-essential gfortran git patch pkg-config unzip \
  python3 python3-venv python3-dev libopenblas-pthread-dev liblapack-dev
```

### 2. Extract the final ZIP into a fresh source directory

Save `cavity-qed-md-final.zip` in `~/Downloads`, then run:

```bash
unzip "$HOME/Downloads/cavity-qed-md-final.zip" -d "$HOME/cavity-qed-v4"
cd "$HOME/cavity-qed-v4/cavity-qed-md"
sha256sum -c SHA256SUMS
```

Use a new extraction directory rather than overlaying an older source tree. Every checksum should report `OK`. A GitHub clone can start from its repository root instead of extracting the ZIP.

For the website's `cavity-qed-md-source.tar.gz` download, use this extraction alternative instead of the ZIP commands:

```bash
mkdir -p "$HOME/cavity-qed-v4"
tar -xzf "$HOME/Downloads/cavity-qed-md-source.tar.gz" -C "$HOME/cavity-qed-v4"
cd "$HOME/cavity-qed-v4/cavity-qed-md"
sha256sum -c SHA256SUMS
```

### 3. Build the corrected CBO-HF installation

From the extracted folder containing `install.sh`:

```bash
BUILD_JOBS=8 WITH_DFT=0 bash install.sh "$HOME/opt/et-ipi-v4"
```

Use a fresh installation prefix for this consolidated release. An interrupted installation made by this installer can resume in the same prefix. Run the installer as your normal user.

The installer creates a private Python environment, installs i-PI 3.2.0 and Python dependencies, builds libcint, extracts the matching eT source, automatically applies the complete v4 gradient correction, generates the required Fortran files, compiles eT and runs the interface tests. `BUILD_JOBS=8` controls compilation parallelism; MD threading is set separately. `WITH_DFT=0` is the default; `WITH_DFT=1` additionally builds LibXC for ordinary DFT single points and does not enable validated CBO-DFT dynamics.

### 4. Activate and validate the compiled executable

After the installation command completes successfully:

```bash
source "$HOME/opt/et-ipi-v4/env.sh"
cd "$HOME/opt/et-ipi-v4/workflow"
python3 compare_backends.py
python3 validate_forces.py --all-nuclear --output force_validation_v4
```

Source `env.sh` in every new shell used for this installation. The backend comparison tests the direct executable and launcher routes. The force suite checks all nine nuclear components and the photon component at lambda = 0 and 0.01 and q = 0 and 0.5. It must finish with `PASS for tested force components`. Keep the default absolute tolerance of 1e-5 a.u. If validation fails, inspect the named log before starting dynamics. Choose a fresh output-directory name when repeating a validation run.

### 5. Run the first monomer smoke test

After validation passes, from the installed `workflow` directory:

```bash
python3 prepare_run.py --xyz examples/water.xyz \
  --eckart off --ensemble nve --steps 100 --timestep-fs 0.1 \
  --temperature 300 --seed 12345 -o water_smoke.xml

ET_THREADS=4 ET_CLIENTS=1 \
ET_FREQUENCY=0.0187505 ET_COUPLING=0.01 \
IPI_INPUT=water_smoke.xml RUN_TAG=water_smoke SOCKET_NAME=water_smoke \
bash run.sh
```

Outputs are under `$HOME/opt/et-ipi-v4/workflow/runs/water_smoke/`, including `ipi.log`, `et_driver_0.log` and the trajectory/property files. Use a new run tag and socket name for another run. This is a 10 fs software smoke test. Production observables require timestep, basis, SCF and sampling convergence. The example uses laboratory-z polarization and the tight cc-pVDZ template.

Installation/build logs are in `$HOME/opt/et-ipi-v4/logs/install.log`; executable and gradient-source identities are recorded in `installation.json` under the same prefix. The correction's independent molecular checks are documented in [GRADIENT_FIX_V4.md](docs/GRADIENT_FIX_V4.md). The rebuilt eT executable still requires the target-machine validation above.

## Joint nuclear–photon PIMD

The original PIMD example is preserved under `examples/monomer-pimd/`, with its own driver, README, tests, templates and launcher. It has different defaults from the newer workstation workflow: polarization along y, frequency 0.01666 a.u., coupling 0.005 a.u., and aug-cc-pVDZ.

From the repository root, after building the patched eT executable:

```bash
source "$HOME/opt/et-ipi-v4/env.sh"
cd examples/monomer-pimd
python3 -m unittest discover -s tests -v
IPI_INPUT=ipi_cbo_pimd.xml ET_CLIENTS=8 ET_THREADS=1 bash run.sh
```

Both nuclei and the photon coordinate have eight beads in this supplied input. Eight beads are a demonstration setting; establish bead convergence. Its PILE-G trajectory samples equilibrium and is not a physical real-time spectrum. The legacy example uses the fixed socket `h2o_cbohf`; run one copy at a time unless its XML and launcher are edited consistently.

## Documentation

| File | Purpose |
|---|---|
| [Method](docs/METHOD.md) | Coordinates, force conventions and scientific scope |
| [Installation](docs/INSTALLATION.md) | Complete installation and workstation reference |
| [Validation](docs/VALIDATION.md) | Checks performed here and real-engine acceptance steps |
| [Provenance](docs/RELEASE_PROVENANCE.md) | Original archive hashes and unchanged source inventory |
| [Publishing](PUBLISHING.md) | Publish the repository and activate GitHub Pages |
| [Licensing](LICENSES.md) | Preserved upstream license and remaining release decision |

PLUMED coupling is reported as implemented in the project discussion, but its current integration files were not present in the recovered bundles. It is not advertised as an executable example in this snapshot. CBO-DFT, correlated-method extensions, and reaction benchmarks are outside this distribution.

## Cite and credit

Credit the eT, i-PI and integral-library developers alongside the interface and any associated scientific publication. Preserve upstream author notices. No software DOI or preprint identifier has been assigned in these files. Add the agreed author list, release license, and permanent identifier before the formal public release; see [PUBLISHING.md](PUBLISHING.md).
