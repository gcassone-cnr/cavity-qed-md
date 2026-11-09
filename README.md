<p align="center">
  <img src="docs/assets/cavity-qed-md-icon.png" width="160" alt="Cavity QED MD: a water molecule and standing wave between cavity mirrors">
</p>

# Cavity-QED Molecular Dynamics

**An eT–i-PI interface for cavity Born–Oppenheimer molecular dynamics and joint nuclear–photon path-integral sampling.**

[Project website](https://gcassone-cnr.github.io/cavity-qed-md/) · [Installation](START_HERE.md) · [Source download](https://gcassone-cnr.github.io/cavity-qed-md/downloads/cavity-qed-md-source.tar.gz) · [Cite the arXiv article](#cite-and-credit)

Maintained by Giuseppe Cassone.

## Features

- Fixed-displacement CBO-HF with analytic nuclear and photon forces.
- Classical NVT and NVE for water monomers and clusters, with one-bead global Eckart NVE.
- Joint nuclear–photon PIMD with PILE-L thermostatting.
- SCF continuation and recovery, OpenMP in eT, and concurrent bead evaluations.
- Force checks, backend comparisons, and analysis of dipoles, spectra and energy conservation.

## Install and validate

The installer targets Ubuntu 22.04 or 24.04 on an x86-64 CPU with Python 3.10–3.12, GNU compilers and LP64 OpenBLAS. It downloads pinned dependencies and builds the supplied eT source with the v4 gradient correction.

### 1. Install system prerequisites

```bash
sudo apt update
sudo apt install -y build-essential gfortran git patch pkg-config unzip \
  python3 python3-venv python3-dev libopenblas-pthread-dev liblapack-dev
```

### 2. Extract the source archive

Download [cavity-qed-md-source.tar.gz](https://gcassone-cnr.github.io/cavity-qed-md/downloads/cavity-qed-md-source.tar.gz) to `~/Downloads`, then run:

```bash
mkdir -p "$HOME/cavity-qed-v4"
tar -xzf "$HOME/Downloads/cavity-qed-md-source.tar.gz" -C "$HOME/cavity-qed-v4"
cd "$HOME/cavity-qed-v4/cavity-qed-md"
sha256sum -c SHA256SUMS
```

Use a fresh source directory and confirm that every checksum reports `OK`. For a Git checkout, run the build commands from the repository root.

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

The monomer PIMD example is under `examples/monomer-pimd/`, with its own driver, README, tests, templates and launcher. It has different defaults from the newer workstation workflow: polarization along y, frequency 0.01666 a.u., coupling 0.005 a.u., and aug-cc-pVDZ.

From the repository root, after building the patched eT executable:

```bash
source "$HOME/opt/et-ipi-v4/env.sh"
cd examples/monomer-pimd
python3 -m unittest discover -s tests -v
IPI_INPUT=ipi_cbo_pimd.xml ET_CLIENTS=8 ET_THREADS=1 bash run.sh
```

Both nuclei and the photon coordinate have eight beads in this supplied input. It uses **PILE-L** (`mode='pile_l'`, `tau=25 fs`), giving the photon centroid its own local Langevin thermalization. Eight beads are a demonstration setting; establish bead convergence. Its thermostatted trajectory samples equilibrium and is not a physical real-time spectrum. Start a fresh run with the updated XML; an old checkpoint retains its stored thermostat. The legacy example uses the fixed socket `h2o_cbohf`; run one copy at a time unless its XML and launcher are edited consistently.

## Documentation

| Guide | Contents |
|---|---|
| [Installation](docs/INSTALLATION.md) | Dependencies, build and example runs |
| [Method](docs/METHOD.md) | CBO-HF model, coordinates and force conventions |
| [Validation](docs/VALIDATION.md) | Electronic force checks and software regression tests |
| [Licensing](LICENSES.md) | License and attribution information |

The supplied workflows cover finite molecular systems with one explicit cavity mode on a CBO-HF surface.

## Cite and credit

If you use **Cavity-QED Molecular Dynamics** in your work, please cite this arXiv article.

> Giuseppe Cassone. *Cavity-QED Molecular Dynamics: An eT-i-PI Interface for Ab Initio Cavity Born-Oppenheimer Dynamics and Joint Nuclear-Photon Path-Integral Sampling.* [arXiv:2611.09667](https://arxiv.org/abs/2611.09667) (2026).

```bibtex
@misc{cassone2026cavityqedmd,
  author = {Cassone, Giuseppe},
  title = {{Cavity-QED Molecular Dynamics}: An {eT}-{i-PI} Interface for Ab Initio Cavity {Born-Oppenheimer} Dynamics and Joint {Nuclear-Photon} Path-Integral Sampling},
  year = {2026},
  eprint = {2611.09667},
  archivePrefix = {arXiv},
  primaryClass = {physics.chem-ph},
  url = {https://arxiv.org/abs/2611.09667}
}
```

[Download BibTeX](docs/cavity-qed-md.bib). GitHub's **Cite this repository** panel uses the arXiv article as the preferred citation. [Citation metadata](CITATION.cff).

Please also cite the relevant eT, i-PI and integral-library publications and preserve upstream author notices. Citation metadata does not change the distribution's licensing; see [LICENSES.md](LICENSES.md).
