#!/usr/bin/env python3
"""Regression check: a clean eT snapshot must have all CMake sources generated.

Uses the exact generation block from install.sh and the bundled eT snapshot.
Runs it twice to cover resuming an interrupted installation; no compiler needed.
"""
from pathlib import Path
import hashlib
import re
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    installer = (ROOT/'install.sh').read_text()
    marker = 'python3 - "$prefix/src/eT-development" <<\'PY\''
    code = installer[installer.index(marker):].split('\n', 1)[1].split('\nPY\n', 1)[0]
    assert installer.index(marker) < installer.index('cmake "${cmake_args[@]}"')
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(ROOT/'vendor/eT-development.zip') as z:
            for name in z.namelist():
                if not (Path(tmp)/name).resolve().is_relative_to(Path(tmp).resolve()):
                    raise ValueError(name)
            z.extractall(tmp)
        source = Path(tmp)/'eT-development'
        missing_file = source/'src/eri_tool/generated_eri_tool/eri_tool_c_class.F90'
        assert not missing_file.exists(), 'Regression fixture already contains generated code'
        subprocess.run(['bash', str(ROOT/'workflow/apply_et_patch.sh'), str(source)], check=True)
        patched = source/'src/wavefunctions/qed_hf/qed_hf_class.F90'
        patched_hash = hashlib.sha256(patched.read_bytes()).hexdigest()
        # Model the already-installed v3 source, then upgrade it back to v4.
        upgrade = ROOT/'workflow/patches/et-fixed-coherent-gradient-v3-to-v4.patch'
        subprocess.run(['patch', '--fuzz=0', '--reverse', '-p1', '-d', str(source)],
                       input=upgrade.read_text(), text=True, check=True)
        subprocess.run(['bash', str(ROOT/'workflow/apply_et_patch.sh'), str(source)], check=True)
        assert hashlib.sha256(patched.read_bytes()).hexdigest() == patched_hash
        subprocess.run(['bash', str(ROOT/'workflow/apply_et_patch.sh'), str(source)], check=True)
        assert hashlib.sha256(patched.read_bytes()).hexdigest() == patched_hash
        for attempt in (1, 2):
            subprocess.run([sys.executable, '-', str(source)], input=code, text=True, check=True)
            inputs = (source/'cmake/eT_sources.cmake',
                      source/'cmake/integrals/libcint/eT_libcint_sources.cmake')
            paths = set()
            for p in inputs:
                paths.update(re.findall(r'(?m)^\s+(src/[^\s)]+\.F90)', p.read_text()))
            assert paths, 'No CMake source paths checked'
            missing = [p for p in sorted(paths) if not (source/p).is_file()]
            assert not missing, f'CMake sources missing: {missing}'
            assert hashlib.sha256(patched.read_bytes()).hexdigest() == patched_hash
            print(f'Generation attempt {attempt}: {len(paths)} CMake sources exist; QED-HF patch unchanged')


if __name__ == '__main__':
    main()
