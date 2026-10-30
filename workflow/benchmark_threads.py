#!/usr/bin/env python3
"""Compare repeated warm eT evaluations on the same short sequence of nearby geometries."""
import argparse
import json
import os
from pathlib import Path
import statistics
import time
import numpy as np
from et_cbohf_ipi_driver import build_parser, ETEvaluator, parse_xyz, flatten, BOHR_PER_ANGSTROM
from launch import physical_cores, ROOT


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xyz',type=Path,default=ROOT/'examples/water21.xyz')
    p.add_argument('--template',type=Path,default=ROOT/'templates/et_cbohf_template_ccpvdz_tight.inp')
    p.add_argument('--threads',type=int,nargs='+')
    p.add_argument('--warm-calls',type=int,default=3)
    p.add_argument('--frequency',type=float,default=0.0187505)
    p.add_argument('--coupling',type=float,default=0.01)
    p.add_argument('--output',type=Path,default=Path('thread_benchmark'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    if a.warm_calls<2:p.error('Use at least two warm calls')
    cores=physical_cores()
    choices=a.threads or sorted(set([n for n in [1,2,4,8,16,32,64] if n<=cores]+[cores]))
    if any(n<1 or n>cores for n in choices):p.error(f'Thread counts must be between 1 and {cores}')
    for k,v in {'OMP_DYNAMIC':'FALSE','OMP_PLACES':'cores','OMP_PROC_BIND':'spread','OMP_MAX_ACTIVE_LEVELS':'1',
                'OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','BLIS_NUM_THREADS':'1'}.items():os.environ[k]=v
    names,rows=parse_xyz(a.xyz);base=np.array([v*BOHR_PER_ANGSTROM for v in flatten(rows[:-1])]+[0.5,0,0])
    results=[];reference=None
    for n in choices:
        args=build_parser().parse_args(['--et-binary',os.environ['ET_BINARY'],'--basis-dir',os.environ['ET_BASIS_DIR'],
              '--template',str(a.template.resolve()),'--threads',str(n),'--frequency',str(a.frequency),
              '--coupling',str(a.coupling),'--polarization','0','0','1','--workdir',str(a.output/f'threads_{n}')])
        evaluator=ETEvaluator(args,names[:-1]);times=[];records=[]
        print(f'Testing {n} threads...',flush=True)
        for i in range(a.warm_calls+1):
            xyz=base.copy();xyz[0]+=i*1e-4
            started=time.perf_counter();energy,forces,metadata=evaluator.evaluate(xyz.tolist());wall=time.perf_counter()-started
            if i:times.append(wall)
            records.append({'energy':energy,'forces':forces,'metadata':metadata,'wall_s':wall})
        if reference is None:reference=records
        error=max(float(np.max(np.abs(np.array(r['forces'])-np.array(s['forces'])))) for r,s in zip(records,reference))
        result={'threads':n,'median_warm_seconds':statistics.median(times),'max_force_difference_vs_first_threads_au':error,'calls':records}
        results.append(result)
        (a.output/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        print(f'{n} threads: {result["median_warm_seconds"]:.3f} s/warm call; max force difference {error:.3g} a.u.',flush=True)
    best=min(results,key=lambda x:x['median_warm_seconds'])
    print(f'Fastest measured: ET_THREADS={best["threads"]}. These nearby-geometry timings are a benchmark, not an MD trajectory.')


if __name__=='__main__':main()
