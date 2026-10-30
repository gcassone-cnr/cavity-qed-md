# CBO-HF / i-PI Eckart recovery bundle

Recovery date: 2026-08-21

This archive collects the latest saved versions of the CBO-HF/i-PI Eckart and
peak-position utilities produced in the associated ChatGPT conversation.

## Provenance and reconstruction status

| File | Status | Notes |
| --- | --- | --- |
| `prepare_eckart_restart.py` | Recovered verbatim | Exact latest saved conversation artifact. |
| `i_pi_eckart_compat.py` | Recovered verbatim | Exact latest saved `.py` artifact; contains both i-PI compatibility fixes described below. |
| `plot_peak_position_diagram.py` | Recovered verbatim | Exact saved molecular/IR peak-position plotting script. |
| `plot_q_peak_position_diagram.py` | Recovered verbatim | Exact saved photon-coordinate peak-position plotting script. |
| `README.md` | Reconstructed | Newly written for this recovery archive from the scripts' own help/docstrings and the validation record in the conversation. |

No Python source file in this archive was recreated from a screenshot or
approximation. The four scripts were recovered from their saved artifacts and
copied without source changes.

## Requirements

- Python 3.10 or later.
- `numpy` for the restart utilities.
- `numpy`, `scipy`, `pandas`, and `matplotlib` for the plotting utilities.
- An activated i-PI environment when using `i_pi_eckart_compat.py`.

For example:

```bash
python3 -m pip install numpy scipy pandas matplotlib
```

## 1. Prepare an aligned Eckart-constrained NVE restart

The default atom ordering is `O H H L`. The script aligns the water bisector
with the laboratory `+z` polarization, places the molecular plane in `xz`,
rotates nuclear momenta consistently, removes the water COM and rigid-body
angular momenta, preserves `q=Lx` and `p_q`, zeros unused `Ly/Lz` components,
and installs an Eckart-constrained NVE motion block.

```bash
python3 prepare_eckart_restart.py RESTART \
  --output ipi_cbo_ir_nve_eckart.xml
```

The input is not overwritten unless `--overwrite` is explicitly supplied.
Inspect the numerical before/after checks printed by the script before launch.

Useful options include:

```bash
python3 prepare_eckart_restart.py --help
```

The generated restart is for one classical bead (`nbeads=1`). It intentionally
uses constrained NVE. The current constrained-dynamics path cannot be combined
with the fixed Cartesian `Ly/Lz` DOFs used in the ordinary input.

## 2. Launch i-PI with the constrained-dynamics compatibility layer

Make the wrapper executable:

```bash
chmod +x i_pi_eckart_compat.py
```

Then use it as the i-PI command in the existing project launcher:

```bash
source "$HOME/venvs/ipi/bin/activate"

IPI_COMMAND="$PWD/i_pi_eckart_compat.py" \
IPI_INPUT=ipi_cbo_ir_nve_eckart.xml \
ET_CLIENTS=1 \
ET_THREADS=8 \
ET_FREQUENCY=0.0187505 \
ET_COUPLING=0.01 \
ET_LAUNCH=/home/cassone/eT-development/build-gcc/eT_launch.py \
RUN_TAG="eckart_NVE_$(date +%Y%m%d_%H%M%S)" \
./run.sh
```

The wrapper applies two compatibility corrections in memory and does not edit
the i-PI installation:

1. the `fixatoms` / `fixatoms_dof` constructor mismatch;
2. the obsolete `forces_mts(level)` call, replaced by the current
   `mts_forces[level].f` interface.

With affected i-PI releases, startup should report both corrections:

```text
i-pi-eckart-compat: applied constrained-dynamics compatibility: fixed-DOF parameter alias, MTS force accessor
```

This wrapper was previously validated in the conversation with a complete
five-step Eckart-constrained NVE test using i-PI 3.2.0 and an internal harmonic
force provider. It propagated and exited normally.

## 3. Molecular/IR peak-position diagram

Each `--spectrum` argument is `OMEGA_C_AU=FILE`. Input may be comma-,
semicolon-, or whitespace-separated; the first two numeric columns are used by
default.

```bash
python3 plot_peak_position_diagram.py \
  --spectrum 0.0186582=Freq0.0186582/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187038=Freq0.0187038/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187505=Freq0.0187505/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0187949=Freq0.0187949/cavity_IR_lab_z_spectrum.csv \
  --spectrum 0.0188404=Freq0.0188404/cavity_IR_lab_z_spectrum.csv \
  --window 4050 4270 \
  --bare-mode 4115 \
  --prominence 0.03 \
  --distance-cm 12 \
  --max-peaks 4 \
  --output water_peak_position.png
```

The output is a two-panel figure containing the detected-maxima validation
panel and the peak-position diagram. A peak table is written by default as
`water_peak_position_peaks.csv`.

## 4. Photon-coordinate q peak-position diagram

### Directly from i-PI trajectories

Repeat the same omega value for independent replicas. `--dt-fs` is the time
between saved XYZ frames, not necessarily the MD integration timestep.

```bash
python3 plot_q_peak_position_diagram.py \
  --trajectory 0.0186582=Freq0.0186582/Sim1/positions_au.xyz \
  --trajectory 0.0186582=Freq0.0186582/Sim2/positions_au.xyz \
  --trajectory 0.0186582=Freq0.0186582/Sim3/positions_au.xyz \
  --dt-fs 0.25 \
  --max-correlation-fs 3000 \
  --broadening-cm 10 \
  --window 4060 4160 \
  --bare-mode 4115.6 \
  --prominence 0.02 \
  --output q_peak_position.png
```

The default photon label is `L`. With `--q-component auto`, the Cartesian
coordinate with the largest variance is selected independently for each file.
Repeated replicas are averaged at the autocorrelation-function level before
the spectrum and maxima are calculated.

### From precomputed q spectra

```bash
python3 plot_q_peak_position_diagram.py \
  --spectrum 0.0186582=Freq0.0186582/q_spectrum.csv \
  --spectrum 0.0187038=Freq0.0187038/q_spectrum.csv \
  --spectrum 0.0187505=Freq0.0187505/q_spectrum.csv \
  --spectrum 0.0187949=Freq0.0187949/q_spectrum.csv \
  --spectrum 0.0188404=Freq0.0188404/q_spectrum.csv \
  --window 4060 4160 \
  --bare-mode 4115.6 \
  --prominence 0.02 \
  --output q_peak_position.png
```

The default trajectory spectrum is the cosine transform of the normalized
`q` autocorrelation. `--observable qdot` changes relative spectral weights but
not the expected peak positions.

## Verification and checksums

All four recovered Python files pass `python3 -m py_compile` in the recovery
environment. The plotting scripts also complete their `--help` paths.

```text
5bf65ae6172633d00fc2f95037d0ce7aee33fd8cbdbe4db25604ef02decf5998  prepare_eckart_restart.py
5f773b7014815dff49575c8b58523249d0ff109a9ae39d3df6eb22c0294180ee  i_pi_eckart_compat.py
2288e617e0fe363c7b6a0fc3c8b62462ad1c31e530b4e50dd87e40fae6de1337  plot_peak_position_diagram.py
ad47340806c81f236ad79fee3f3490ce85bcab317e8a34f10bbe7b300bffcfe1  plot_q_peak_position_diagram.py
```

