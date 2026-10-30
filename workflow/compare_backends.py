#!/usr/bin/env python3
"""Check direct execution against the original eT launcher on the same geometry."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from et_cbohf_ipi_driver import ETEvaluator, build_parser, parse_xyz, flatten, BOHR_PER_ANGSTROM
from launch import ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xyz',type=Path,default=ROOT/'examples/water.xyz')
    p.add_argument('--template',type=Path,default=ROOT/'templates/et_cbohf_template_ccpvdz_tight.inp')
    p.add_argument('--frequency',type=float,default=0.0187505)
    p.add_argument('--coupling',type=float,default=0.01)
    p.add_argument('--q',type=float,default=0.5)
    p.add_argument('--threads',type=int,default=1)
    p.add_argument('--output',type=Path,default=Path('backend_validation'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    names,rows=parse_xyz(a.xyz)
    if names[-1]!='L' or names.count('L')!=1:p.error('Need exactly one final L')
    positions=[v*BOHR_PER_ANGSTROM for v in flatten(rows[:-1])]+[a.q,0,0]
    records={}
    for backend in ('launcher','direct'):
        argv=['--et-launch',os.environ['ET_LAUNCH'],'--template',str(a.template),
              '--frequency',str(a.frequency),'--coupling',str(a.coupling),'--polarization','0','0','1',
              '--threads',str(a.threads),'--workdir',str(a.output/backend)]
        if backend=='direct':argv+=['--et-binary',os.environ['ET_BINARY'],'--basis-dir',os.environ['ET_BASIS_DIR']]
        ev=ETEvaluator(build_parser().parse_args(argv),names[:-1])
        records[backend]=[]
        for label in ('cold','warm'):
            energy,forces,metadata=ev.evaluate(positions)
            records[backend].append({'label':label,'energy':energy,'forces':forces,'metadata':metadata})
    (a.output/'results.json').write_text(json.dumps(records,indent=2)+'\n')
    for left,right in zip(records['launcher'],records['direct']):
        de=abs(left['energy']-right['energy'])
        df=float(np.max(np.abs(np.array(left['forces'])-right['forces'])))
        dm=float(np.max(np.abs(np.array(left['metadata']['dipole_au'])-right['metadata']['dipole_au'])))
        print(f'{left["label"]}: |dE|={de:.3g} Eh; max |dF|={df:.3g} au; max |dmu|={dm:.3g} au')
        if de>1e-9 or df>1e-8 or dm>1e-8:
            raise SystemExit('FAIL: investigate backend differences before using direct execution')
    print('PASS: direct/launcher agreement on the tested cold and warm evaluations')


if __name__=='__main__':main()
