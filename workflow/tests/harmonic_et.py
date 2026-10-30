#!/usr/bin/env python3
"""TEST DOUBLE ONLY: analytic harmonic potential with a coupled photon; not eT."""
import math
import os
from pathlib import Path
import re
import sys

t = Path('eT.inp').read_text()
if (Path.cwd().parent/'always_fail').exists():
    sys.exit(10)
restart=Path('orbital_coefficients_qed-rhf')
if restart.exists() and restart.read_text().strip()=='corrupt':
    sys.exit(9)
def scalar(key):
    return float(re.search(key+r':\s*\{([^}]+)\}',t).group(1))
w=scalar('frequency');lam=scalar('coupling');q=scalar('coherent state')/math.sqrt(w/2)
pol=[float(v) for v in re.search(r'polarization:\s*\{([^}]+)\}',t).group(1).split(',')]
atoms=re.findall(r'^\s*(O|H)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)',t,re.M)
mu=[sum((-2 if row[0]=='O' else 1)*float(row[i+1]) for row in atoms) for i in range(3)]
mup=sum(x*y for x,y in zip(mu,pol))
k=0.0001
e=-76+sum(k*float(v)**2 for a in atoms for v in a[1:])+0.5*(w*q+lam*mup)**2
iterations=2 if restart.exists() else 13
text=f'Total energy: {e:.16e}\nConvergence criterion met in {iterations} iterations!\nDipole moment [Debye]:\nComp. Electronic Nuclear Total\n'
for axis,m in zip('xyz',mu):
    text+=f'{axis} 0.0000000000 0.0000000000 {m/0.393430270:.16e}\n'
text+='eT terminated successfully!\n'
Path('eT.out').write_text(text)
if 'molecular gradient' in t:
    out=' - Molecular gradient (Hartree/bohr)\n Center x y z\n'
    for i,a in enumerate(atoms):
        c=-2 if a[0]=='O' else 1
        g=[2*k*float(a[j+1])+(w*q+lam*mup)*lam*c*pol[j] for j in range(3)]
        out+=f'{i+1} {a[0]} '+' '.join(f'{v:.16e}' for v in g)+'\n'
    Path('eT.molecular_gradient.out').write_text(out)
restart.write_text('valid\n')
Path('observed_environment.txt').write_text('\n'.join(k+'='+os.environ.get(k,'') for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','AO_BASIS_PATH']))
