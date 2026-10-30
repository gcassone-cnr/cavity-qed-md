# Complete corrected source package — v4

This folder contains the matching eT source archive, the eT–i-PI drivers,
MD/PIMD examples, analysis tools, documentation and all corrections discussed
in the installation session. No earlier package or separately downloaded
repair patch is needed.

The installer automatically applies the complete v4 fixed-displacement
nuclear-gradient patch and generates eT's required Fortran sources before
CMake configuration. The default installation is CBO-HF.

This is a source distribution for Ubuntu 22.04/24.04 and Python 3.10–3.12.
Compilers and system libraries are required; Internet access is needed to
fetch the pinned Python dependencies and libcint. It is not a precompiled
or offline installation image.

## Installation steps

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

## Included verification

- 27 interface tests passed; one full socket test is skipped where UNIX
  sockets are prohibited.
- Clean source generation, upgrade from v3 and repeated v4 application passed.
- An independent molecular CBO-HF calculation reproduced the reported old
  nuclear-force mismatch; the corrected algebra agreed with energy differences
  below 1e-9 a.u. for nonzero-coupling q = 0 and 0.5 monomer checks.

The maintainer confirmed the corrected v4 installation and documented test workflow
as working on 30 October 2026. Repeat the target-machine checks in step 4 for each
installation. See `docs/GRADIENT_FIX_V4.md` and `docs/VALIDATION.md` for evidence
and limits; raw workstation logs are not bundled.
