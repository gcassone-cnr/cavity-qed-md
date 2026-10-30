> Adapted from the original 18 September installation document for the prepared source tree. The original text is retained as `ORIGINAL_INSTALLATION.md` at the repository root.

**Corrected v4 installation — 30 October 2026**

The README and `START_HERE.md` at the repository root provide the numbered installation and first-run steps. The sections below add workstation and cluster details.

This bundle assembles the recovered eT source, its fixed-coherent-state gradient patch, the optimized CBO-HF driver, and the recovered i-PI Eckart compatibility wrapper. The new launcher supports a monomer or a 21-water cluster without hard-coded atom indices. It provides classical, one-bead global Eckart NVE and ordinary NVT/NVE workflows.

The installer targets **Ubuntu 22.04/24.04 on an x86-64 CPU workstation**, with Python 3.10–3.12, GNU compilers, and OpenBLAS. It needs Internet access to download the pinned dependencies. It includes the matching eT source archive; no old workstation files are required. It builds for the local CPU, so rebuild on a different workstation rather than copying its executable.

**Validation status:** workflow regression tests and real i-PI 3.2.0 harmonic Eckart tests were run here. The electronic eT build and real electronic forces could not be tested in this execution environment because a Fortran compiler was unavailable. UNIX socket creation was also prohibited; the complete socket smoke test is included and runs during installation on the workstation. Read `validation/REPORT.md` for the exact scope. Installation success alone is not a validation of CBO forces or production trajectories.

**1. Install the system prerequisites**

On the new workstation:

```bash
sudo apt update
sudo apt install -y build-essential gfortran git patch pkg-config unzip \
  python3 python3-venv python3-dev libopenblas-pthread-dev liblapack-dev
```

Use the LP64 OpenBLAS package, not an ILP64/64-bit-integer replacement. The installer sets eT's default integers to 32 bits to match it. This follows eT's [build configuration requirements](https://etprogram.org/setup.html).

**2. Build from the prepared source tree**

For the final ZIP saved in `~/Downloads`, extract into a fresh source directory:

```bash
unzip "$HOME/Downloads/cavity-qed-md-final.zip" -d "$HOME/cavity-qed-v4"
cd "$HOME/cavity-qed-v4/cavity-qed-md"
```

For the website's tar.gz download, use this alternative:

```bash
mkdir -p "$HOME/cavity-qed-v4"
tar -xzf "$HOME/Downloads/cavity-qed-md-source.tar.gz" -C "$HOME/cavity-qed-v4"
cd "$HOME/cavity-qed-v4/cavity-qed-md"
```

From the chosen extracted source directory, verify and build:

```bash
sha256sum -c SHA256SUMS
BUILD_JOBS=8 WITH_DFT=0 bash install.sh "$HOME/opt/et-ipi-v4"
source "$HOME/opt/et-ipi-v4/env.sh"
cd "$HOME/opt/et-ipi-v4/workflow"
```

`BUILD_JOBS` controls compilation parallelism, not MD threads. Reduce it if compilation runs out of memory. Use a fresh installation prefix. Re-running after a failed build can resume a prefix created by this installer; it does not overwrite an unrelated installation. Use a new prefix for a different bundle or compiler.

The installer builds libcint 6.1.2 and eT in Release mode with OpenMP, using the native CPU instruction set. It installs i-PI 3.2.0 in a private virtual environment. It applies the complete v4 eT patch automatically and generates the required Fortran sources before CMake configuration. There is no separate manual patch step. The default is `WITH_DFT=0`; setting `WITH_DFT=1` also builds LibXC 5.1.3 for ordinary restricted DFT single-point calculations.

The generated `env.sh` records the executable, basis directory and i-PI paths. Source it in every new shell. Installation logs, linked libraries and Python versions are retained under the installation's `logs/`; `installation.json` records source identities.

The corrected installer runs eT's own source generator before configuring CMake. This produces the complex Fortran modules and submodule interfaces absent from the source archive. A clean-snapshot and repeat-generation regression check is available with `python3 scripts/check_et_generation.py` from the source-tree root. If an earlier installer stopped at a missing `generated_eri_tool/eri_tool_c_class.F90`, use the corrected package and rerun the installation command with the same managed installation prefix; the existing dependencies can be reused.

**3. Validate the installed executable before dynamics**

```bash
python3 -m unittest discover -s tests -v
python3 compare_backends.py
python3 validate_forces.py --all-nuclear --output force_validation_v4
```

The backend comparison checks both cold and restarted calculations through the original launcher and the faster direct-binary route. It uses the monomer, so this should be much cheaper than debugging a cluster trajectory.

The force validation checks all nine monomer nuclear components and the photon force at q=0 and 0.5, with coupling 0 and 0.01. The default central-difference step is 1e-4 bohr for nuclei and 1e-4 atomic units for q; the absolute force tolerance is 1e-5 atomic units. These are numerical checks, not claims that every parameter regime has been verified. If a check fails, inspect the logged errors, repeat with a smaller finite-difference step, and diagnose the discrepancy before MD. Do not increase the tolerance simply to obtain a pass.

For the actual 21-water configuration and chosen cavity parameters, also run:

```bash
python3 validate_forces.py --xyz examples/water21.xyz \
  --threads 16 --frequency 0.0187505 --couplings 0.01 \
  --output force_validation_21
```

Choose a thread count available on your workstation. Without `--all-nuclear`, this checks three nuclear components and the photon component at both q values. For broader validation, add `--all-nuclear`, use representative trajectory configurations, and include q values spanning the trajectory with `--q`. Full checks on 63 nuclei are expensive. Use the same template, frequency and coupling intended for production.

**4. Run a short 21-water Eckart test**

The recovered `examples/water21.xyz` contains 63 real nuclei and a final `L` particle. To use your own configuration, supply an XYZ with 64 rows in the desired nuclear order and `L` last. Nuclear coordinates are in angstrom; **L.x is q in atomic units**, while L.y and L.z are zero. The photon mass is set to one electron mass. The molecular orientation is retained; cavity polarization is laboratory z.

Generate an explicit input, then launch it:

```bash
python3 prepare_run.py --xyz examples/water21.xyz \
  --eckart global --ensemble nve --steps 20 --timestep-fs 0.1 \
  --temperature 300 --seed 12345 -o water21_smoke.xml

ET_THREADS=auto ET_CLIENTS=1 \
ET_TEMPLATE=templates/et_cbohf_template_ccpvdz_tight.inp \
ET_FREQUENCY=0.0187505 ET_COUPLING=0.01 \
IPI_INPUT=water21_smoke.xml RUN_TAG=water21_smoke \
SOCKET_NAME=water21_smoke ./run.sh
```

Frequency and coupling above are explicit example parameters, not an optimized or universally appropriate cavity. Replace them consistently in force validation, benchmarking and MD. The 20-step run is a software smoke test, not an IR production trajectory. Fresh random thermal momenta do not make an equilibrated cluster.

`./run.sh --dry-run` stages and checks the input without running either engine. Every actual launch first checks one nuclear and the photon force at q=0 and 0.5. For 21 waters, the photon index is 63 (zero-based). Unique run tags and socket names are required for concurrent replicas. A used run tag is never overwritten.

**Eckart behavior:** all 63 nuclei form one global constraint group. Six external degrees of freedom are removed: center-of-mass translation and rotation relative to a fixed reference. Bond lengths, angles and relative molecular motions remain active. The photon is excluded. The reference is explicitly stored and preserved across Eckart checkpoints. This is constrained nuclear dynamics; in a cavity with a fixed polarization it changes the allowed rotational dynamics and should be used as part of the specified physical model.

For ordinary NVT/NVE, L.y/L.z are fixed as Cartesian DOFs 190 and 191. For Eckart NVE, those Cartesian constraints are omitted because the recovered i-PI constrained path does not support that combination. Their coordinates and momenta start at zero and remain zero because the driver returns zero transverse photon forces. Tests check this behavior. Keep `fixcom=false` with Eckart. This bundle supports Eckart **one-bead NVE**; it does not advertise Eckart NVT or Eckart PIMD as validated features.

**5. Equilibrate and prepare production segments**

To prepare ordinary NVT at 300 K:

```bash
python3 prepare_run.py --xyz examples/water21.xyz \
  --eckart off --ensemble nvt --steps 2000 --timestep-fs 0.1 \
  --temperature 300 -o water21_nvt.xml
IPI_INPUT=water21_nvt.xml RUN_TAG=water21_nvt SOCKET_NAME=water21_nvt ./run.sh
```

The 2000 steps are an example duration; establish equilibration from the trajectory. An unconstrained NVT ensemble and a constrained one have different degrees of freedom. Switching to Eckart removes external kinetic energy, and the preparation script reports how much. Discard the resulting transient rather than assuming that the projected state is already a stationary constrained ensemble.

To start a new Eckart NVE segment from its final i-PI checkpoint:

```bash
python3 prepare_run.py --restart runs/water21_nvt/RESTART \
  --eckart global --ensemble nve --steps 10000 --timestep-fs 0.1 \
  -o water21_nve.xml
IPI_INPUT=water21_nve.xml RUN_TAG=water21_nve SOCKET_NAME=water21_nve ./run.sh
```

Use the actual final `RESTART` path, or the most recent complete `traj.checkpoint`, after inspecting the run directory. Conversion preserves nuclear positions and photon q/p_q, projects nuclear momenta into the allowed subspace, resets the step counter and regenerates output/run-control fields. If the source was already a global Eckart checkpoint, its original reference is retained. This creates a new segment; it is not a bit-for-bit continuation of every integrator or random-number-generator field.

For a direct continuation of an existing input/checkpoint, set `IPI_INPUT` to that file and increase its `total_steps` beyond the stored step. The new run uses a new output directory and may start with cold electronic orbitals. The launcher does not silently reset a checkpoint's step counter.

**6. Select the fastest measured thread count**

```bash
python3 benchmark_threads.py --xyz examples/water21.xyz \
  --threads 8 16 32 --warm-calls 3 --output thread_benchmark
```

Only request counts available on the machine; add 64 on a 64-physical-core workstation. Each count performs a cold call followed by three nearby-geometry warm calls, with identical settings. The script saves energies, forces, SCF metadata and wall times and prints the fastest measured count. Use it as `ET_THREADS=N` for subsequent runs. The benchmark is not an MD trajectory or a guarantee of production throughput. Inspect force differences across counts as well as timing.

`ET_THREADS=auto` uses physical cores within the process's CPU affinity and divides them among clients. More threads can be slower because of memory bandwidth and NUMA effects. For one classical bead use `ET_CLIENTS=1`; extra clients do not advance a serial trajectory concurrently. Independent replicas can run in parallel with disjoint CPU allocations and unique sockets.

Measure production throughput and energy behavior with:

```bash
python3 utilities/summarize_driver_performance.py \
  runs/water21_nve/et_driver_0.log --timestep-fs 0.1 --discard 2
python3 utilities/summarize_nve.py \
  runs/water21_nve/traj.properties --waters 21
```

The first utility estimates simulated fs/day from force-evaluation times; startup and other overhead can reduce actual throughput. The second reports conserved-energy drift and fluctuations. Compare matched NVE segments at 0.1 and 0.05 fs. Consider 0.2 fs only after assessing drift and spectral changes at equal physical duration. No timestep increase is automatically approved by this bundle.

**7. Performance settings and scientific choices**

| Setting | Default | Effect |
|---|---|---|
| `ET_BACKEND` | `direct` | Runs the native eT executable without an extra Python launcher or launcher scratch-copy cycle. |
| `ET_THREADS` | `auto` | Physical-core allocation; override using measured timings. |
| `ET_BLAS_THREADS` | `1` | Controls OpenBLAS/MKL/BLIS threads within the OpenMP calculation. |
| `ET_SCF_RESTART` | `auto` | Uses the previous converged orbitals; retries cold if the warm call fails. |
| `ET_KEEP_EVERY` | `0` | Retains only the current calculation, valid orbitals and failures. Set e.g. 1000 for sparse snapshots. |
| `ET_SCRATCH_ROOT` | Run directory | Set to a local SSD directory to reduce network filesystem traffic. |
| `ET_TEMPLATE` | Tight cc-pVDZ | Basis/SCF model; select explicitly when comparing with older runs. |
| `ET_FORCE_TOL` | `1e-5` | Absolute tolerance for mandatory launch force checks, in atomic units. |

For example, `ET_SCRATCH_ROOT="$HOME/et-scratch"` places electronic work under that directory, separated by run tag. i-PI trajectories and checkpoints stay under `workflow/runs/`. Avoid a network-mounted home for electronic scratch when a local SSD is available. Do not delete scratch during a run.

The driver retains `current/`, `scf_restart/`, optional `snapshots/`, and failure diagnostics for each client. Restart candidates are promoted only after successful output validation; failed warm/cold calls do not replace the last valid orbitals. Startup finite-difference displacements use energy-only inputs. The electronic executable still starts for each force call and rebuilds integrals; this is not a persistent in-process electronic engine.

All four templates are included:

| Template under `templates/` | Basis | Energy / residual thresholds |
|---|---|---|
| `et_cbohf_template.inp` | aug-cc-pVDZ | 1e-10 / 1e-9 |
| `et_cbohf_template_aug_relaxed.inp` | aug-cc-pVDZ | 1e-8 / 1e-7 |
| `et_cbohf_template_ccpvdz_tight.inp` | cc-pVDZ | 1e-10 / 1e-9 |
| `et_cbohf_template_production_ccpvdz.inp` | cc-pVDZ | 1e-8 / 1e-7 |

The default uses cc-pVDZ with tight SCF thresholds. **Changing from aug-cc-pVDZ to cc-pVDZ changes the electronic model.** To reproduce an earlier aug-cc-pVDZ run, set `ET_TEMPLATE=templates/et_cbohf_template.inp`. Basis and SCF changes require snapshot force/dipole comparisons, fresh equilibration as appropriate, and NVE validation. Threading, direct execution and scratch retention do not intentionally change the Hamiltonian.

`utilities/compare_probe_results.py` compares driver `--probe-json` results, including force RMSE, maximum force error and dipole error. Its optional numerical gates are user-selected accuracy targets. Passing a few snapshots does not establish ensemble or spectral equivalence.

**8. DFT scope**

With `WITH_DFT=1`, test ordinary PBE with:

```bash
mkdir -p dft_single_point
cp examples/pbe_single_point.inp dft_single_point/
cd dft_single_point
"$ET_LAUNCH" pbe_single_point.inp --omp 8
cd ..
```

This is an ordinary restricted DFT energy calculation, not cavity DFT MD. The [official eT DFT documentation](https://etprogram.org/dft_calculation.html) describes restricted DFT single points. Source inspection and recovery did not establish a validated fixed-q CBO-DFT nuclear-gradient implementation. The CBO driver therefore rejects a DFT template. Replacing `qed-hf` with `dft`, or adding LibXC, does not supply the missing cavity energy/force model.

A future CBO-DFT extension must define its energy functional and self-polarization terms, implement matching nuclear and photon derivatives, and pass finite-difference, lambda=0, restart and NVE tests. Ordinary DFT speed does not by itself establish correct cavity forces. No CBO-DFT acceleration is claimed for this bundle.

**9. Analysis and troubleshooting**

The recovered `analysis/` directory includes IR, photon-coordinate spectra, replica averaging and peak plots. Use `--help` on the individual tools. Output spacing is `timestep × trajectory stride`; the generated inputs use stride 1 unless changed. For example:

```bash
python3 analysis/ir_spectrum.py runs/water21_nve/traj.dipole_z_au_0 \
  --sample-dt-fs 0.1 --output-prefix runs/water21_nve/ir_z
```

Use the actual extras filename emitted by i-PI, and discard equilibration appropriately. The legacy monomer-specific Eckart preparation script is retained under `legacy_eckart/` for provenance; the current cluster workflow uses `prepare_run.py`.

- Patch mismatch: use the bundled eT snapshot. Do not force the patch onto unrelated source.
- LibXC build failure: inspect `logs/install.log`; use `WITH_DFT=0` if you only need the HF route.
- Missing i-PI or wrong Python: source the generated `env.sh` again. The pip distribution is named **ipi**, pinned to 3.2.0.
- Warm SCF failure: inspect the per-client `failures/` snapshot; the driver automatically retries from a cold guess. Repeated fallback can erase the restart speed benefit.
- Force-check failure: inspect `force_check_q0.log` or `force_check_q05.log`. A source patch reporting “already applied” does not prove the executable was rebuilt or its forces are correct.
- Socket already exists: choose a new socket name and check for a still-running job. The supervisor only stops process groups that it started.
- Existing run directory: choose a new `RUN_TAG`. The script does not overwrite prior trajectories.
- NVE drift: compare a smaller timestep at the same physical duration, check SCF convergence and force errors, and inspect constraint residuals. A short smoke test cannot certify long production dynamics.

Source identities and changes are recorded in `PROVENANCE.md`. The bundled eT source retains its GPLv3 license and author notices. Dependency licensing remains with each upstream project.
