#!/usr/bin/env python3
"""Stage, force-check and supervise an eT/i-PI run. Settings use environment variables."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import numpy as np
from prepare_run import array, names_array, BOHR_PER_ANGSTROM

ROOT = Path(__file__).resolve().parent


def physical_cores():
    cpus = os.sched_getaffinity(0)
    cores = set()
    for cpu in cpus:
        topology = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        try:
            cores.add(((topology/'physical_package_id').read_text().strip(), (topology/'core_id').read_text().strip()))
        except OSError:
            cores.add(('cpu', str(cpu)))
    return max(len(cores), 1)


def resolve(value):
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (ROOT/p).resolve()


def stage(source, dest, socket, clients):
    tree = ET.parse(source)
    root = tree.getroot()
    systems = root.findall('system')
    if len(systems) != 1:
        raise ValueError('This launcher supports exactly one system')
    system = systems[0]
    beads = system.find('beads')
    if beads is None:
        raise ValueError('Use prepare_run.py to create an input with embedded coordinates/momenta')
    names = names_array(beads.findtext('names', ''))
    nbeads = int(beads.get('nbeads', '0'))
    if nbeads < 1 or int(beads.get('natoms', '0')) != len(names):
        raise ValueError('Bead counts and atom names are inconsistent')
    for tag in ('q','p','m'):
        value = beads.find(tag)
        if value is None or value.get('units','atomic_unit') != 'atomic_unit':
            raise ValueError('Embedded q/p/m must be in atomic units; use prepare_run.py')
    if names.count('L') != 1 or names[-1] != 'L':
        raise ValueError('Expected exactly one final L photon')
    q = array(beads.findtext('q', '')).reshape(nbeads, len(names), 3)
    p = array(beads.findtext('p', '')).reshape(nbeads, len(names), 3)
    m = array(beads.findtext('m', ''))
    if len(m) != len(names) or np.any(m <= 0) or not np.isclose(m[-1], 1.0):
        raise ValueError('Photon mass must be 1 in atomic units')
    if any(not np.isfinite(x).all() for x in (q,p,m)):
        raise ValueError('Nonfinite input coordinates, momenta or masses')
    if np.any(q[:, -1, 1:] != 0) or np.any(p[:, -1, 1:] != 0):
        raise ValueError('L.y/L.z coordinates and momenta must start at zero')
    motion = system.find('motion')
    constrained = motion is not None and motion.get('mode') == 'constrained_dynamics'
    if constrained:
        cd = motion.find('constrained_dynamics')
        if nbeads != 1 or cd is None or cd.get('mode') != 'nve':
            raise ValueError('Recovered Eckart route supports one-bead NVE only')
        fixed = motion.find('fixatoms_dof')
        if fixed is not None and array(fixed.text or '').size:
            raise ValueError('Eckart input must not contain fixed Cartesian DOFs')
        if motion.findtext('fixcom', 'false').strip().lower() in ('true','1'):
            raise ValueError('Eckart requires fixcom=false')
        constraints = cd.findall('constraint')
        expected = list(range(len(names)-1))
        if len(constraints) != 1 or constraints[0].get('mode') != 'eckart' or array(constraints[0].findtext('atoms', '')).astype(int).tolist() != expected:
            raise ValueError('This bundle supports one global Eckart group of all real nuclei')
    else:
        if motion is None or motion.get('mode') != 'dynamics':
            raise ValueError('Expected dynamics or constrained_dynamics')
        if motion.findtext('fixcom', 'false').strip().lower() in ('true','1'):
            raise ValueError('Use fixcom=false: the generic COM includes the photon pseudo-atom')
        fixed = motion.find('fixatoms_dof')
        if fixed is None:
            fixed = ET.SubElement(motion, 'fixatoms_dof')
        indices = set(array(fixed.text or '').astype(int).tolist())
        indices.update([3*len(names)-2, 3*len(names)-1])
        fixed.text = '[ ' + ', '.join(map(str, sorted(indices))) + ' ]'
    ff = root.findall('ffsocket')
    if len(ff) != 1 or ff[0].get('mode') != 'unix':
        raise ValueError('Expected one UNIX ffsocket')
    ff[0].set('pbc', 'false')
    ff[0].find('address').text = socket
    slots = ff[0].find('slots')
    if slots is None:
        slots = ET.SubElement(ff[0], 'slots')
    slots.text = str(clients)
    tree.write(dest/'input.xml', encoding='unicode')
    with (dest/'probe.xyz').open('w') as f:
        f.write(f'{len(names)}\nNuclei: angstrom; L.x: photon q in atomic units\n')
        for s, r in zip(names[:-1], q[0, :-1]/BOHR_PER_ANGSTROM):
            f.write(s+' '+' '.join(f'{v:.16e}' for v in r)+'\n')
        f.write('L '+' '.join(f'{v:.16e}' for v in q[0,-1])+'\n')
    return names, constrained


def stop(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run', action='store_true', help='Stage and validate the input without running eT/i-PI')
    a = parser.parse_args()
    clients = int(os.environ.get('ET_CLIENTS', '1'))
    if clients < 1:
        raise ValueError('ET_CLIENTS must be positive')
    setting = os.environ.get('ET_THREADS', 'auto')
    threads = max(1, physical_cores()//clients) if setting == 'auto' else int(setting)
    if threads < 1:
        raise ValueError('ET_THREADS must be positive')
    tag = os.environ.get('RUN_TAG', time.strftime('%Y%m%d_%H%M%S'))
    socket = os.environ.get('SOCKET_NAME', 'et_'+tag)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', tag) or not re.fullmatch(r'[A-Za-z0-9_.-]{1,80}',socket):
        raise ValueError('Run tag/socket must use only letters, digits, dot, underscore or hyphen')
    run = resolve(os.environ.get('RUN_ROOT', 'runs')) / tag
    run.mkdir(parents=True, exist_ok=False)
    source = resolve(os.environ.get('IPI_INPUT', 'examples/water_eckart.xml'))
    names, constrained = stage(source, run, socket, clients)
    environment = os.environ.copy()
    for key, value in {'OMP_DYNAMIC':'FALSE','OMP_PLACES':'cores','OMP_PROC_BIND':'spread','OMP_MAX_ACTIVE_LEVELS':'1'}.items():
        environment.setdefault(key, value)
    for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
        environment[key] = os.environ.get('ET_BLAS_THREADS', '1')
    template = resolve(os.environ.get('ET_TEMPLATE', 'templates/et_cbohf_template_ccpvdz_tight.inp'))
    backend = os.environ.get('ET_BACKEND', 'direct')
    common = [sys.executable, str(ROOT/'et_cbohf_ipi_driver.py'), '--template', str(template),
              '--threads', str(threads), '--frequency', os.environ.get('ET_FREQUENCY','0.0187505'),
              '--coupling', os.environ.get('ET_COUPLING','0.01'), '--polarization','0','0','1',
              '--scf-restart', os.environ.get('ET_SCF_RESTART','auto'),
              '--keep-every', os.environ.get('ET_KEEP_EVERY','0'), '--nuclear-symbols', *names[:-1]]
    if backend == 'direct':
        if not a.dry_run and (not os.environ.get('ET_BINARY') or not os.environ.get('ET_BASIS_DIR')):
            raise ValueError('Direct execution needs ET_BINARY and ET_BASIS_DIR; source the installed env.sh')
        common += ['--et-binary', os.environ.get('ET_BINARY', ''), '--basis-dir', os.environ.get('ET_BASIS_DIR','')]
    elif backend == 'launcher':
        common += ['--et-launch', os.environ.get('ET_LAUNCH','eT_launch.py')]
    else:
        raise ValueError('ET_BACKEND must be direct or launcher')
    if not template.is_file():
        raise FileNotFoundError(template)
    (run/'launch.json').write_text(json.dumps({'source':str(source),'nuclei':len(names)-1,'eckart':constrained,
                   'clients':clients,'threads_per_client':threads,'command':common,'environment':{k:v for k,v in environment.items() if k.startswith(('OMP_','OPENBLAS_','MKL_','BLIS_'))}}, indent=2)+'\n')
    print(f'Run: {run}\nNuclei: {len(names)-1}; global Eckart: {constrained}; clients: {clients}; threads/client: {threads}',flush=True)
    if a.dry_run:
        return
    socket_path = Path('/tmp/ipi_'+socket)
    if socket_path.exists():
        raise FileExistsError(f'Socket already exists: {socket_path}; choose another SOCKET_NAME')
    scratch = resolve(os.environ.get('ET_SCRATCH_ROOT', str(run/'et_runs'))) 
    if 'ET_SCRATCH_ROOT' in os.environ:
        scratch = scratch / tag
    scratch.mkdir(parents=True, exist_ok=False)
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    children, handles = [], []
    try:
        for label, q in [('q0','0'),('q05','0.5')]:
            print(f'Checking forces at q={q} ...',flush=True)
            log = (run/f'force_check_{label}.log').open('w'); handles.append(log)
            command = common + ['--probe',str(run/'probe.xyz'),'--probe-q',q,
                       '--check-force','0,z','--check-force',f'{len(names)-1},x',
                       '--fd-tolerance',os.environ.get('ET_FORCE_TOL','1e-5'),
                       '--probe-json',str(run/f'force_check_{label}.json'),
                       '--workdir',str(scratch/f'check_{label}')]
            process = subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
            children.append(process)
            if process.wait() != 0:
                raise RuntimeError(f'Force check failed: {log.name}')
        ipi = [sys.executable, str(ROOT/'i_pi_eckart_compat.py')] if constrained else [os.environ.get('IPI_COMMAND','i-pi')]
        log = (run/'ipi.log').open('w'); handles.append(log)
        server = subprocess.Popen(ipi+['input.xml'], cwd=run,stdout=log,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
        children.append(server)
        deadline = time.monotonic()+60
        while not socket_path.is_socket():
            if server.poll() is not None:
                raise RuntimeError(f'i-PI exited before opening its socket; inspect {log.name}')
            if time.monotonic()>deadline:
                raise RuntimeError(f'Timeout opening {socket_path}; inspect {log.name}')
            time.sleep(0.1)
        drivers=[]
        for i in range(clients):
            log=(run/f'et_driver_{i}.log').open('w');handles.append(log)
            process=subprocess.Popen(common+['--unix',socket,'--workdir',str(scratch/f'client_{i}')],cwd=run,stdout=log,stderr=subprocess.STDOUT,env=environment,start_new_session=True)
            children.append(process);drivers.append(process)
        while server.poll() is None:
            for i,process in enumerate(drivers):
                if process.poll() is not None:
                    # i-PI closes clients during finalization, before exiting.
                    if process.returncode == 0:
                        try:
                            server.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            pass
                    if server.poll() is None:
                        raise RuntimeError(f'Client {i} exited before i-PI; inspect et_driver_{i}.log')
            time.sleep(0.1)
        if server.returncode != 0:
            raise RuntimeError(f'i-PI failed; inspect {run}/ipi.log')
        for process in drivers:
            try:
                rc=process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                raise RuntimeError('i-PI finished but an eT client did not exit')
            if rc != 0:
                raise RuntimeError('eT client reported an error; inspect its log')
        print(f'Completed: {run}',flush=True)
    finally:
        for process in reversed(children):
            if process.poll() is None:
                stop(process)
        for log in handles:
            log.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as error:
        print(f'ERROR: {error}',file=sys.stderr)
        sys.exit(1)
