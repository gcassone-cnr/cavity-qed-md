#!/usr/bin/env bash
set -euo pipefail

project_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"

: "${ET_LAUNCH:?Set ET_LAUNCH to the path of eT_launch.py from the patched build}"
IPI_COMMAND="${IPI_COMMAND:-i-pi}"
IPI_INPUT="${IPI_INPUT:-ipi_cbo_nve.xml}"
ET_THREADS="${ET_THREADS:-1}"
ET_CLIENTS="${ET_CLIENTS:-1}"
ET_FREQUENCY="${ET_FREQUENCY:-0.01666}"
ET_COUPLING="${ET_COUPLING:-0.005}"
ET_GRADIENT_MODE="${ET_GRADIENT_MODE:-analytic}"
run_tag="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
run_dir="${project_dir}/runs/${run_tag}"

if ! command -v "${IPI_COMMAND}" >/dev/null 2>&1; then
  echo "Cannot find i-PI command: ${IPI_COMMAND}" >&2
  exit 1
fi
if ! command -v "${ET_LAUNCH}" >/dev/null 2>&1 && [[ ! -x "${ET_LAUNCH}" ]]; then
  echo "Cannot execute eT launcher: ${ET_LAUNCH}" >&2
  exit 1
fi
if [[ ! -f "${project_dir}/${IPI_INPUT}" ]]; then
  echo "Cannot find i-PI input: ${project_dir}/${IPI_INPUT}" >&2
  exit 1
fi
if [[ ! "${ET_CLIENTS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ET_CLIENTS must be a positive integer" >&2
  exit 1
fi
if [[ -e "${run_dir}" ]]; then
  echo "Run directory already exists: ${run_dir}" >&2
  exit 1
fi

mkdir -p "${run_dir}"
cp "${project_dir}/${IPI_INPUT}" "${run_dir}/input.xml"
cp "${project_dir}/h2o_cbo_init.xyz" "${run_dir}/"

common_args=(
  --et-launch "${ET_LAUNCH}"
  --template "${project_dir}/et_cbohf_template.inp"
  --threads "${ET_THREADS}"
  --frequency "${ET_FREQUENCY}"
  --coupling "${ET_COUPLING}"
  --gradient-mode "${ET_GRADIENT_MODE}"
)

echo "Run directory: ${run_dir}"
echo "Checking nuclear and photon forces at q=0 against energy differences..."
python3 "${project_dir}/et_cbohf_ipi_driver.py" \
  "${common_args[@]}" \
  --probe "${run_dir}/h2o_cbo_init.xyz" \
  --check-force 0,y \
  --check-force 3,x \
  --workdir "${run_dir}/et_force_check_q0" \
  >"${run_dir}/force_check_q0.log" 2>&1

echo "Checking nuclear and photon forces at q=0.5 au against energy differences..."
python3 "${project_dir}/et_cbohf_ipi_driver.py" \
  "${common_args[@]}" \
  --probe "${run_dir}/h2o_cbo_init.xyz" \
  --probe-q 0.5 \
  --check-force 0,y \
  --check-force 3,x \
  --workdir "${run_dir}/et_force_check_q05" \
  >"${run_dir}/force_check_q05.log" 2>&1

cd "${run_dir}"
"${IPI_COMMAND}" input.xml >ipi.log 2>&1 &
ipi_pid=$!
driver_pids=()

cleanup() {
  for pid in "${driver_pids[@]:-}"; do
    if kill -0 "${pid}" 2>/dev/null; then
      kill "${pid}" 2>/dev/null || true
    fi
  done
  if kill -0 "${ipi_pid}" 2>/dev/null; then
    kill "${ipi_pid}" 2>/dev/null || true
  fi
}
trap cleanup EXIT INT TERM

socket_path="/tmp/ipi_h2o_cbohf"
socket_ready=false
for _ in {1..300}; do
  if ! kill -0 "${ipi_pid}" 2>/dev/null; then
    wait "${ipi_pid}" || true
    echo "i-PI exited before opening ${socket_path}. Last log lines:" >&2
    tail -n 40 ipi.log >&2
    exit 1
  fi
  if [[ -S "${socket_path}" ]]; then
    socket_ready=true
    break
  fi
  sleep 0.1
done
if [[ "${socket_ready}" != true ]]; then
  echo "Timed out waiting for i-PI socket ${socket_path}" >&2
  exit 1
fi

for ((client = 0; client < ET_CLIENTS; client++)); do
  python3 "${project_dir}/et_cbohf_ipi_driver.py" \
    "${common_args[@]}" \
    --unix h2o_cbohf \
    --workdir "${run_dir}/et_runs/client_${client}" \
    >"et_driver_${client}.log" 2>&1 &
  driver_pids+=("$!")
done

wait "${ipi_pid}"
for pid in "${driver_pids[@]}"; do
  wait "${pid}"
done
trap - EXIT INT TERM
echo "Completed: ${run_dir}"
