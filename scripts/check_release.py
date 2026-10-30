#!/usr/bin/env python3
"""Check source identity, the downloadable archive, and local website links."""
from pathlib import Path
import ast
import hashlib
import subprocess
import tarfile
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]

def check_manifest(path):
    count = 0
    for line in path.read_text().splitlines():
        digest, rel = line.split('  ', 1)
        p = path.parent/rel
        if hashlib.sha256(p.read_bytes()).hexdigest() != digest:
            raise AssertionError(f'Checksum mismatch: {p.relative_to(ROOT)}')
        count += 1
    return count

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links = []; self.ids = set()
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if 'id' in a:
            if a['id'] in self.ids: raise AssertionError(f'Duplicate id {a["id"]}')
            self.ids.add(a['id'])
        for key in ('href', 'src'):
            if key in a: self.links.append(a[key])

def main():
    count = check_manifest(ROOT/'SHA256SUMS')
    check_manifest(ROOT/'docs/downloads/SHA256SUMS')
    with tarfile.open(ROOT/'docs/downloads/cavity-qed-md-source.tar.gz') as t:
        members = {m.name:m for m in t.getmembers()}
        for line in (ROOT/'SHA256SUMS').read_text().splitlines():
            digest, rel = line.split('  ', 1)
            member = members['cavity-qed-md/'+rel]
            if not member.isfile(): raise AssertionError(member.name)
            if hashlib.sha256(t.extractfile(member).read()).hexdigest() != digest:
                raise AssertionError(f'Stale source download: {rel}')
    page = ROOT/'docs/index.html'
    links = Links(); links.feed(page.read_text())
    for url in links.links:
        s = urlsplit(url)
        if s.scheme or s.netloc: continue
        if not s.path:
            if s.fragment and s.fragment not in links.ids: raise AssertionError(url)
        elif not (page.parent/unquote(s.path)).is_file():
            raise AssertionError(f'Missing page resource: {url}')
    for p in ROOT.rglob('*.py'):
        if '__pycache__' not in p.parts: ast.parse(p.read_text(), filename=str(p))
    for p in ROOT.rglob('*.sh'):
        subprocess.run(['bash', '-n', str(p)], check=True)
    print(f'Passed: {count} source checksums, source download contents, page links/anchors, Python and Bash syntax.')

if __name__ == '__main__':
    main()
