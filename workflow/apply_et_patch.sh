#!/usr/bin/env bash
set -euo pipefail

project_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
et_source="${1:-${ET_SOURCE:-}}"
patch_file="${project_dir}/patches/et-fixed-coherent-gradient.patch"
upgrade_v1_v2="${project_dir}/patches/et-fixed-coherent-gradient-v1-to-v2.patch"
upgrade_v2_v3="${project_dir}/patches/et-fixed-coherent-gradient-v2-to-v3.patch"
upgrade_v3_v4="${project_dir}/patches/et-fixed-coherent-gradient-v3-to-v4.patch"
target="src/wavefunctions/qed_hf/qed_hf_class.F90"

if [[ -z "${et_source}" ]]; then
  echo "Usage: $0 /absolute/path/to/eT-source" >&2
  exit 2
fi
if [[ ! -f "${et_source}/${target}" ]]; then
  echo "Not an eT source tree (missing ${target}): ${et_source}" >&2
  exit 2
fi

if patch --fuzz=0 --dry-run --forward -p1 -d "${et_source}" <"${patch_file}" >/dev/null; then
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${patch_file}"
  echo "Applied fixed-coherent-state gradient patch to ${et_source}"
elif patch --fuzz=0 --dry-run --reverse -p1 -d "${et_source}" <"${patch_file}" >/dev/null; then
  echo "Current fixed-coherent-state gradient patch is already applied in ${et_source}"
elif patch --fuzz=0 --dry-run --forward -p1 -d "${et_source}" <"${upgrade_v3_v4}" >/dev/null; then
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v3_v4}"
  echo "Upgraded fixed-coherent-state gradient patch from v3 to v4 in ${et_source}"
elif patch --fuzz=0 --dry-run --forward -p1 -d "${et_source}" <"${upgrade_v2_v3}" >/dev/null; then
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v2_v3}"
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v3_v4}"
  echo "Upgraded fixed-coherent-state gradient patch from v2 to v4 in ${et_source}"
elif patch --fuzz=0 --dry-run --forward -p1 -d "${et_source}" <"${upgrade_v1_v2}" >/dev/null; then
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v1_v2}"
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v2_v3}"
  patch --fuzz=0 --forward -p1 -d "${et_source}" <"${upgrade_v3_v4}"
  echo "Upgraded fixed-coherent-state gradient patch from v1 to v4 in ${et_source}"
else
  echo "Patch does not apply cleanly. The eT source may differ from the uploaded development snapshot." >&2
  exit 1
fi
