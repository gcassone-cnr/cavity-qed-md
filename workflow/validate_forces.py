#!/usr/bin/env python3
"""Validate the actual compiled eT nuclear/photon forces before production."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
from et_cbohf_ipi_driver import parse_xyz

ROOT=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xyz',type=Path,default=ROOT/'examples/water.xyz')
    p.add_argument('--template',type=Path,default=ROOT/'templates/et_cbohf_template_ccpvdz_tight.inp')
    p.add_argument('--all-nuclear',action='store_true')
    p.add_argument('--q',type=float,nargs='+',default=[0.0,0.5])
    p.add_argument('--couplings',type=float,nargs='+',default=[0.0,0.01])
    p.add_argument('--frequency',type=float,default=0.0187505)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--tolerance',type=float,default=1e-5)
    p.add_argument('--fd-step',type=float,default=1e-4)
    p.add_argument('--output',type=Path,default=Path('force_validation'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    names,_=parse_xyz(a.xyz)
    if names.count('L')!=1 or names[-1]!='L':
        p.error('Need exactly one final L')
    indices=[f'{i},{axis}' for i in range(len(names)-1) for axis in 'xyz'] if a.all_nuclear else ['0,z','1,x','1,z']
    indices += [f'{len(names)-1},x']
    for li,coupling in enumerate(a.couplings):
        for qi,q in enumerate(a.q):
            tag=f'lambda{li}_q{qi}'
            cmd=[sys.executable,str(ROOT/'et_cbohf_ipi_driver.py'),'--et-binary',os.environ['ET_BINARY'],
                 '--basis-dir',os.environ['ET_BASIS_DIR'],'--template',str(a.template.resolve()),
                 '--probe',str(a.xyz.resolve()),'--nuclear-symbols','auto','--polarization','0','0','1',
                 '--probe-q',str(q),'--frequency',str(a.frequency),'--coupling',str(coupling),
                 '--threads',str(a.threads),'--fd-tolerance',str(a.tolerance),'--fd-check-step',str(a.fd_step),
                 '--workdir',str(a.output/tag),'--probe-json',str(a.output/(tag+'.json'))]
            for index in indices:cmd+=['--check-force',index]
            print(f'Validating lambda={coupling}, q={q}',flush=True)
            with (a.output/(tag+'.log')).open('w') as f:
                result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
            if result.returncode:
                raise SystemExit(f'FAILED: inspect {a.output/(tag+".log")}. Resolve the discrepancy; do not loosen the tolerance just to pass.')
    print(f'PASS for tested force components; results in {a.output}. Next: NVE timestep/drift validation.')


if __name__=='__main__':main()
