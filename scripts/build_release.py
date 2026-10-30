#!/usr/bin/env python3
"""Refresh the checksums and the self-contained website source download."""
from pathlib import Path
import hashlib
import io
import tarfile
import gzip

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {'.git', '__pycache__', 'runs', '.venv', 'venv', 'downloads'}

def source_files():
    for p in sorted(ROOT.rglob('*')):
        rel = p.relative_to(ROOT)
        if p.is_symlink():
            raise ValueError(f'Symlinks are not part of this release: {rel}')
        if not p.is_file() or any(x in EXCLUDED for x in rel.parts):
            continue
        if p.name in {'SHA256SUMS', '.DS_Store'} or p.suffix in {'.pyc', '.pyo'}:
            continue
        if p.name == '.env' or p.name.startswith('.env.'):
            continue
        yield p, rel.as_posix()

def main():
    files = list(source_files())
    manifest = ''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {r}\n' for p, r in files)
    (ROOT/'SHA256SUMS').write_text(manifest)
    dest = ROOT/'docs/downloads'
    dest.mkdir(parents=True, exist_ok=True)
    archive = dest/'cavity-qed-md-source.tar.gz'
    with archive.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as gz, tarfile.open(fileobj=gz, mode='w') as t:
        for p, rel in files + [(ROOT/'SHA256SUMS', 'SHA256SUMS')]:
            data = p.read_bytes()
            info = tarfile.TarInfo('cavity-qed-md/'+rel)
            info.size = len(data)
            info.mode = 0o755 if p.stat().st_mode & 0o111 else 0o644
            info.mtime = 0
            t.addfile(info, io.BytesIO(data))
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (dest/'SHA256SUMS').write_text(f'{digest}  {archive.name}\n')
    print(f'Prepared {len(files)} files; source download {archive.stat().st_size:,} bytes; SHA256 {digest}')

if __name__ == '__main__':
    main()
