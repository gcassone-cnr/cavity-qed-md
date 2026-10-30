#!/usr/bin/env bash
# Ubuntu 22.04 / 24.04, GNU compilers, OpenBLAS LP64.
set -euo pipefail
bundle_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
prefix="${1:-${HOME}/opt/et-ipi-workstation}"
build_jobs="${BUILD_JOBS:-8}"
with_dft="${WITH_DFT:-0}"
if [[ ! "$build_jobs" =~ ^[1-9][0-9]*$ || ! "$with_dft" =~ ^[01]$ ]]; then
  echo 'BUILD_JOBS must be positive; WITH_DFT must be 0 or 1' >&2
  exit 2
fi
for tool in python3 git gcc g++ gfortran patch; do
  command -v "$tool" >/dev/null || { echo "Missing prerequisite: $tool (see README.md)" >&2; exit 1; }
done
python3 - <<'PY'
import sys
if not (3,10) <= sys.version_info[:2] <= (3,12):
    raise SystemExit('Use Python 3.10, 3.11 or 3.12 for this pinned installation')
PY
prefix="$(python3 -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).expanduser().resolve())' "$prefix")"
if [[ -e "$prefix" && ! -f "$prefix/.et-ipi-installer" ]]; then
  echo "Prefix already exists and is not managed by this installer: $prefix" >&2
  echo 'Choose a fresh prefix.' >&2
  exit 1
fi
mkdir -p "$prefix"
touch "$prefix/.et-ipi-installer"
mkdir -p "$prefix/src" "$prefix/logs"
exec > >(tee -a "$prefix/logs/install.log") 2>&1
echo "Installing to $prefix"
python3 -m venv "$prefix/venv"
"$prefix/venv/bin/python" -m pip install -r "$bundle_dir/requirements.txt"
export PATH="$prefix/venv/bin:$PATH"

# Build the integral library with eT's required ordering and cutoff settings.
if [[ ! -d "$prefix/src/libcint/.git" ]]; then
  git clone --depth 1 --branch v6.1.2 https://github.com/sunqm/libcint.git "$prefix/src/libcint"
fi
[[ "$(git -C "$prefix/src/libcint" describe --tags --exact-match)" == v6.1.2 ]]
cmake -S "$prefix/src/libcint" -B "$prefix/src/libcint/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DPYPZPX=1 -DMIN_EXPCUTOFF=0.0 \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_C_COMPILER=gcc \
  -DCMAKE_C_FLAGS_RELEASE='-O3 -DNDEBUG -march=native' \
  -DCMAKE_INSTALL_PREFIX="$prefix/libcint"
cmake --build "$prefix/src/libcint/build" --parallel "$build_jobs"
cmake --install "$prefix/src/libcint/build"

if [[ "$with_dft" == 1 ]]; then
  if [[ ! -d "$prefix/src/libxc/.git" ]]; then
    git clone --depth 1 --branch 5.1.3 https://gitlab.com/libxc/libxc.git "$prefix/src/libxc"
  fi
  [[ "$(git -C "$prefix/src/libxc" describe --tags --exact-match)" == 5.1.3 ]]
  cmake -S "$prefix/src/libxc" -B "$prefix/src/libxc/build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -DENABLE_FORTRAN=ON -DENABLE_FORTRAN03=ON \
    -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_INSTALL_PREFIX="$prefix/libxc"
  cmake --build "$prefix/src/libxc/build" --parallel "$build_jobs"
  cmake --install "$prefix/src/libxc/build"
fi

# Exact original snapshot, not a moving development branch.
python3 - "$bundle_dir/vendor/eT-development.zip" "$prefix/src" <<'PY'
import hashlib, pathlib, sys, zipfile
archive=pathlib.Path(sys.argv[1]); dest=pathlib.Path(sys.argv[2])
expected='45724c57cf4f3c28c889e0749f650d9aaef49e3c99690361525073b458dbe92d'
if hashlib.sha256(archive.read_bytes()).hexdigest()!=expected:
    raise SystemExit('eT source archive checksum mismatch')
if not (dest/'eT-development').exists():
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            target=(dest/name).resolve()
            if not target.is_relative_to(dest.resolve()):
                raise SystemExit('Unsafe archive member')
        z.extractall(dest)
PY
bash "$bundle_dir/workflow/apply_et_patch.sh" "$prefix/src/eT-development"
# eT's setup.py generates complex sources and submodule interfaces before
# invoking CMake. Direct CMake configuration must perform the same step.
# Pass the source root explicitly: the generator lives under dev_tools/.
python3 - "$prefix/src/eT-development" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / 'dev_tools'))
from autogenerate_files import autogenerate
autogenerate(root)
required = root / 'src/eri_tool/generated_eri_tool/eri_tool_c_class.F90'
if not required.is_file():
    raise SystemExit(f'eT source generation did not create {required}')
print('Generated eT complex sources and submodule interfaces')
PY
export LIBCINT_ROOT="$prefix/libcint"
cmake_args=(-S "$prefix/src/eT-development" -B "$prefix/build-et" -G Ninja
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc
  -DCMAKE_CXX_COMPILER=g++ -DENABLE_OMP=ON -DENABLE_64BIT_INTEGERS=OFF
  -DINTEGRAL_LIBRARY=libcint -DBLAS_TYPE=OPENBLAS -DLAPACK_TYPE=SYSTEM_NATIVE
  -DENABLE_RUNTIME_CHECKS=OFF -DENABLE_FORCED_BATCHING=OFF -DENABLE_INITIALIZE_NAN=OFF)
if [[ "$with_dft" == 1 ]]; then
  export LIBXC_ROOT="$prefix/libxc"
  cmake_args+=(-DLIBXC_ROOT="$LIBXC_ROOT" -DCMAKE_DISABLE_FIND_PACKAGE_LibXC=FALSE)
else
  cmake_args+=(-DCMAKE_DISABLE_FIND_PACKAGE_LibXC=TRUE)
fi
cmake "${cmake_args[@]}"
if [[ "$with_dft" == 1 ]]; then
  python3 - "$prefix/build-et/CMakeCache.txt" <<'PY'
import pathlib, re, sys
cache=pathlib.Path(sys.argv[1]).read_text()
match=re.search(r'^Libxc_DIR:PATH=(.+)$',cache,re.M)
if not match or 'NOTFOUND' in match[1] or not pathlib.Path(match[1]).is_dir():
    raise SystemExit('WITH_DFT=1 but LibXC was not found by eT; inspect the configure log')
PY
fi
cmake --build "$prefix/build-et" --parallel "$build_jobs"

# Keep the workflow outside the eT source. Refresh only on first installation.
if [[ ! -d "$prefix/workflow" ]]; then
  cp -a "$bundle_dir/workflow" "$prefix/workflow"
fi
chmod +x "$prefix/build-et/eT_launch.py" "$prefix/workflow/run.sh" "$prefix/workflow/tests/"*.py
python3 - "$prefix" "$with_dft" <<'PY'
import hashlib,json,pathlib,shlex,subprocess,sys
p=pathlib.Path(sys.argv[1]); q=shlex.quote
lines=['# Source this file in bash before running the workflow.',f'source {q(str(p/"venv/bin/activate"))}']
env={'ET_BINARY':str(p/'build-et/eT'),'ET_LAUNCH':str(p/'build-et/eT_launch.py'),
     'ET_BASIS_DIR':str(p/'src/eT-development/ao_basis'),'IPI_COMMAND':str(p/'venv/bin/i-pi'),
     'ET_BACKEND':'direct','ET_THREADS':'auto','ET_BLAS_THREADS':'1',
     'ET_SCF_RESTART':'auto','ET_KEEP_EVERY':'0','OMP_PLACES':'cores','OMP_PROC_BIND':'spread',
     'OMP_DYNAMIC':'FALSE','OMP_MAX_ACTIVE_LEVELS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1'}
lines += ['export '+k+'='+q(v) for k,v in env.items()]
(p/'env.sh').write_text('\n'.join(lines)+'\n')
record={'eT_original_source_sha256':'45724c57cf4f3c28c889e0749f650d9aaef49e3c99690361525073b458dbe92d',
        'gradient_patch_version':4,
        'eT_patched_gradient_sha256':hashlib.sha256((p/'src/eT-development/src/wavefunctions/qed_hf/qed_hf_class.F90').read_bytes()).hexdigest(),
        'eT_sources_generated':True,
        'ipi':'3.2.0','libcint':subprocess.check_output(['git','-C',str(p/'src/libcint'),'rev-parse','HEAD'],text=True).strip(),
        'libxc_enabled':sys.argv[2]=='1'}
if record['libxc_enabled']:
 record['libxc']=subprocess.check_output(['git','-C',str(p/'src/libxc'),'rev-parse','HEAD'],text=True).strip()
(p/'installation.json').write_text(json.dumps(record,indent=2)+'\n')
PY
python3 -m pip freeze > "$prefix/logs/python-freeze.txt"
ldd "$prefix/build-et/eT" > "$prefix/logs/et-linked-libraries.txt"
cd "$prefix/workflow"
python3 -m unittest discover -s tests -v
echo "Build and workflow tests completed. Next: source '$prefix/env.sh'"
echo "Then cd '$prefix/workflow' and run: python3 validate_forces.py --all-nuclear"
echo 'Real electronic force and NVE validation is still required before production.'
