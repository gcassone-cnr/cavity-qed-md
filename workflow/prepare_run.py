#!/usr/bin/env python3
"""Create a one-bead eT/i-PI input from an angstrom XYZ or an i-PI RESTART.

The last XYZ row is L: its x value is photon q in atomic units, not angstrom.
Global Eckart fixes six external degrees of freedom of all real nuclei.
It does not freeze individual water orientations, bond lengths or angles.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
from ipi.utils.units import Constants, Elements, unit_to_internal

BOHR_PER_ANGSTROM = 1.8897261254578281


def array(text):
    return np.asarray([float(x.replace('D', 'E').replace('d', 'e')) for x in
                       re.findall(r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eEdD][-+]?\d+)?', text)])


def names_array(text):
    return [x.strip("'\"") for x in re.split(r'[\s,\[\]]+', text.strip()) if x]


def fmt(values):
    return '[ ' + ', '.join(f'{x:.16e}' for x in np.asarray(values).ravel()) + ' ]'


def add(parent, tag, text=None, **attrs):
    node = ET.SubElement(parent, tag, attrs)
    if text is not None:
        node.text = str(text)
    return node


def project_external(q, p, m):
    p = p.copy()
    p -= m[:, None] * p.sum(axis=0) / m.sum()
    r = q - np.average(q, axis=0, weights=m)
    inertia = sum(mi * (ri @ ri * np.eye(3) - np.outer(ri, ri)) for mi, ri in zip(m, r))
    if np.linalg.matrix_rank(inertia) != 3:
        raise ValueError('Global Eckart requires a nonlinear group of nuclei')
    omega = np.linalg.solve(inertia, np.cross(r, p).sum(axis=0))
    p -= m[:, None] * np.cross(omega, r)
    return p


def load_xyz(path, temperature, seed):
    from et_cbohf_ipi_driver import parse_xyz
    names, rows = parse_xyz(path)
    if names[-1] != 'L' or names.count('L') != 1:
        raise ValueError('XYZ must contain exactly one final L photon particle')
    q = np.asarray(rows, dtype=float)
    q[:-1] *= BOHR_PER_ANGSTROM
    m = np.asarray([Elements.mass(s) for s in names])
    # A photon mass of one electron mass is essential for the q convention.
    m[-1] = 1.0
    thermal_energy = Constants.kb * unit_to_internal('temperature', 'kelvin', temperature)
    p = np.random.default_rng(seed).normal(size=q.shape) * np.sqrt(m[:, None] * thermal_energy)
    return names, q, p, m


def load_restart(path):
    root = ET.parse(path).getroot()
    b = root.findall('./system/beads')
    if len(b) != 1 or int(b[0].get('nbeads', '0')) != 1:
        raise ValueError('Expected one system and one classical bead')
    b = b[0]
    for tag in ('q', 'p', 'm'):
        if b.find(tag) is None or b.find(tag).get('units', 'atomic_unit') != 'atomic_unit':
            raise ValueError('RESTART q/p/m must be present in atomic units')
    names = names_array(b.findtext('names', ''))
    q = array(b.findtext('q')).reshape(-1, 3)
    p = array(b.findtext('p')).reshape(-1, 3)
    m = array(b.findtext('m'))
    if q.shape != p.shape or len(names) != len(m) or len(q) != len(m):
        raise ValueError('Inconsistent RESTART atom arrays')
    if names[-1] != 'L' or names.count('L') != 1 or not np.isclose(m[-1], 1.0):
        raise ValueError('Expected one final L with mass 1 in atomic units')
    return names, q, p, m


def restart_reference(path, names):
    """Preserve the original Eckart frame when converting an Eckart checkpoint."""
    root = ET.parse(path).getroot()
    motion = root.find('./system/motion')
    if motion is None or motion.get('mode') != 'constrained_dynamics':
        return None
    constraints = motion.findall('./constrained_dynamics/constraint')
    if len(constraints) != 1 or constraints[0].get('mode') != 'eckart':
        raise ValueError('Expected one global Eckart constraint in the checkpoint')
    c = constraints[0]
    if array(c.findtext('atoms', '')).astype(int).tolist() != list(range(len(names)-1)):
        raise ValueError('Checkpoint Eckart group must contain all real nuclei')
    value = c.find('values')
    if value is None or value.get('units', 'atomic_unit') != 'atomic_unit':
        raise ValueError('Checkpoint must store its original Eckart reference in atomic units')
    ref = array(value.text or '')
    if ref.size != 3*(len(names)-1) or not np.isfinite(ref).all():
        raise ValueError('Missing or invalid Eckart reference; use the original input to recover it')
    return ref.reshape(-1, 3)


def build(names, q, p, m, *, eckart, ensemble, steps, timestep, temperature, stride, seed, reference=None):
    if not all(np.isfinite(v).all() for v in (q, p, m)) or np.any(m <= 0):
        raise ValueError('Coordinates, momenta and positive masses must be finite')
    if eckart and ensemble != 'nve':
        raise ValueError('This recovered Eckart route supports NVE only')
    q, p = q.copy(), p.copy()
    q[-1, 1:] = p[-1, 1:] = 0.0
    removed = 0.0
    if eckart:
        reference = q[:-1].copy() if reference is None else np.asarray(reference).copy()
        if reference.shape != q[:-1].shape or not np.isfinite(reference).all():
            raise ValueError('Invalid Eckart reference')
        rref = reference - np.average(reference, axis=0, weights=m[:-1])
        delta = q[:-1] - reference
        residual = np.r_[np.average(delta, axis=0, weights=m[:-1]),
                         np.sum(np.cross(m[:-1,None]*rref, delta), axis=0)/m[:-1].sum()]
        if np.linalg.norm(residual) > 1e-8:
            raise ValueError('Checkpoint coordinates do not satisfy their original Eckart constraints')
        before = np.sum(p[:-1]**2/m[:-1, None])/2
        # The Eckart tangent condition uses the fixed reference, not the
        # instantaneous angular momentum of a deformed configuration.
        p[:-1] = project_external(reference, p[:-1], m[:-1])
        removed = before - np.sum(p[:-1]**2/m[:-1, None])/2
    root = ET.Element('simulation', verbosity='low')
    output = add(root, 'output', prefix='traj')
    add(output, 'properties', '[ step, time{femtosecond}, conserved{electronvolt}, potential{electronvolt}, kinetic_md{electronvolt} ]', stride='1', filename='properties')
    add(output, 'trajectory', 'positions{atomic_unit}', stride=str(stride), filename='positions_au', format='xyz', cell_units='atomic_unit')
    add(output, 'trajectory', 'extras', stride=str(stride), filename='dipole_xyz_au', extra_type='dipole_au')
    add(output, 'trajectory', 'extras', stride=str(stride), filename='dipole_z_au', extra_type='polarized_dipole_au')
    add(output, 'checkpoint', filename='checkpoint', stride='100', overwrite='True')
    add(root, 'total_steps', steps)
    add(add(root, 'prng'), 'seed', seed)
    ff = add(root, 'ffsocket', name='et-cbohf', mode='unix', matching='auto', pbc='false')
    add(ff, 'address', 'et_cbohf')
    add(ff, 'latency', '0.001')
    add(ff, 'timeout', '0')
    add(ff, 'slots', '1')
    system = add(root, 'system')
    b = add(system, 'beads', natoms=str(len(names)), nbeads='1')
    add(b, 'q', fmt(q), shape=f'(1, {3*len(names)})')
    add(b, 'p', fmt(p), shape=f'(1, {3*len(names)})')
    add(b, 'm', fmt(m), shape=f'({len(names)})')
    add(b, 'names', '[ ' + ', '.join(names) + ' ]', shape=f'({len(names)})')
    # Nonperiodic engine: cell is bookkeeping, not an electronic PBC box.
    box = np.max(np.ptp(q[:-1], axis=0)) + 40*BOHR_PER_ANGSTROM
    add(system, 'cell', fmt(np.eye(3)*box), shape='(3, 3)')
    add(add(system, 'forces'), 'force', forcefield='et-cbohf')
    motion = add(system, 'motion', mode='constrained_dynamics' if eckart else 'dynamics')
    add(motion, 'fixcom', 'false')
    if not eckart:
        add(motion, 'fixatoms_dof', f'[ {3*len(names)-2}, {3*len(names)-1} ]')
    dynamics = add(motion, 'constrained_dynamics' if eckart else 'dynamics', mode=ensemble)
    add(dynamics, 'timestep', timestep, units='femtosecond')
    if ensemble == 'nvt':
        add(add(dynamics, 'thermostat', mode='langevin'), 'tau', '100', units='femtosecond')
    if eckart:
        add(dynamics, 'nsteps_geo', '5')
        cs = add(dynamics, 'csolver')
        add(cs, 'tolerance', '1e-10')
        add(cs, 'maxit', '1000')
        add(cs, 'norm_order', '2')
        constraint = add(dynamics, 'constraint', mode='eckart')
        add(constraint, 'atoms', '[ ' + ', '.join(map(str, range(len(names)-1))) + ' ]')
        # Store a fixed reference explicitly, including across RESTARTs.
        add(constraint, 'values', fmt(reference), units='atomic_unit')
    add(add(system, 'ensemble'), 'temperature', temperature, units='kelvin')
    return root, removed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--xyz', type=Path)
    source.add_argument('--restart', type=Path)
    parser.add_argument('-o', '--output', required=True, type=Path)
    parser.add_argument('--eckart', choices=['global', 'off'], default='global')
    parser.add_argument('--ensemble', choices=['nve', 'nvt'], default='nve')
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--timestep-fs', type=float, default=0.1)
    parser.add_argument('--temperature', type=float, default=300)
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--seed', type=int, default=314159)
    a = parser.parse_args()
    if a.output.exists():
        parser.error('Output exists; choose a new output name')
    if min(a.steps, a.timestep_fs, a.temperature, a.stride) <= 0:
        parser.error('Steps, timestep, temperature, and stride must be positive')
    names, q, p, m = load_xyz(a.xyz, a.temperature, a.seed) if a.xyz else load_restart(a.restart)
    reference = restart_reference(a.restart, names) if a.restart and a.eckart=='global' else None
    root, removed = build(names, q, p, m, eckart=a.eckart=='global', ensemble=a.ensemble, steps=a.steps,
                          timestep=a.timestep_fs, temperature=a.temperature, stride=a.stride, seed=a.seed,
                          reference=reference)
    ET.indent(root)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(a.output, encoding='unicode')
    print(f'Wrote {a.output}: {len(names)-1} nuclei + L; {a.eckart} Eckart; {a.ensemble}')
    print(f'Removed external kinetic energy: {removed:.12g} hartree; preserved photon q and p_q')
    if a.xyz:
        print('Fresh thermal momenta are not an equilibrated ensemble. Discard initial equilibration.')
    else:
        print('Starts a new segment at step zero; retains the stored Eckart reference when present.')


if __name__ == '__main__':
    main()
