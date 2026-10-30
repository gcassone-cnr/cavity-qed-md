"""Optional independent CBO-HF gradient check using PySCF 2.14.0.

Install PySCF in a separate test environment; it is not an eT runtime dependency.
The one-electron DSE uses the same finite-basis mu S^-1 mu convention as the
bundled eT. Nuclear gradients combine PySCF's ordinary integral/Pulay derivative
with cavity terms; S and mu derivatives use central differences. Full relaxed
energy differences provide the reference. This checks the fixed-displacement
algebra and reproduces the reported water failure, not compiled eT gradients.
"""
from pathlib import Path
import numpy as np
from pyscf import gto, scf, lib

lib.num_threads(1)
ANGSTROM_TO_BOHR = 1.8897261254578281
ROOT = Path(__file__).resolve().parents[1]
rows = (ROOT/'workflow/examples/water.xyz').read_text().splitlines()[2:5]
symbols = [r.split()[0] for r in rows]
coords = np.array([[float(v) for v in r.split()[1:]] for r in rows])*ANGSTROM_TO_BOHR
omega = .0187505
lam = .01

def mol_at(r):
    m = gto.M(atom=list(zip(symbols, r)), unit='Bohr', basis='cc-pvdz', verbose=0)
    m.incore_anyway = True  # This small reference fits in memory.
    return m

def matrices(m):
    s = m.intor_symmetric('int1e_ovlp')
    mu_nuc = np.dot(m.atom_charges(), m.atom_coords()[:,2])
    mu = -m.intor_symmetric('int1e_r', comp=3)[2] + mu_nuc/m.nelectron*s
    return s, mu

def solve(r, q, coupling=lam):
    m = mol_at(r)
    s,mu = matrices(m)
    mf = scf.RHF(m)
    h = mf.get_hcore()
    ordinary_veff = mf.get_veff
    extra_h = .5*coupling**2*mu@np.linalg.solve(s,mu) + omega*coupling*q*mu + .5*omega**2*q*q/m.nelectron*s
    mf.get_hcore = lambda mol=None: h+extra_h
    def veff(mol=None, dm=None, *args, **kwargs):
        if dm is None: dm=mf.make_rdm1()
        return ordinary_veff(mol,dm) + coupling**2*(np.einsum('ij,ji',dm,mu)*mu-.5*mu@dm@mu)
    mf.get_veff = veff
    mf.conv_tol=1e-13
    mf.conv_tol_grad=1e-10
    mf.max_cycle=200
    mf.kernel()
    if not mf.converged: raise RuntimeError('Reference SCF failed')
    return mf

if __name__ == '__main__':
    for coupling,q in [(0.,0.), (0.,.5), (.01,0.), (.01,.5)]:
        mf=solve(coords,q,coupling)
        print('lambda=',coupling,'q=',q,'energy=',mf.e_tot,flush=True)
        s,mu=matrices(mf.mol)
        P=mf.make_rdm1();D=P/2
        si=np.linalg.inv(s);d=coupling*mu
        dexp=np.einsum('ij,ji',P,d)
        A=D@d@(si-D)
        pbase=-si@d@D@d@si
        z=q*np.sqrt(omega/2)
        base=mf.nuc_grad_method().kernel()
        gnew=base.copy();gold=base.copy();gfd=np.zeros((3,3))
        h=1e-4
        for a in range(3):
            for k in range(3):
                plus=coords.copy();minus=coords.copy();plus[a,k]+=h;minus[a,k]-=h
                sp,mp=matrices(mol_at(plus));sm,mm=matrices(mol_at(minus))
                ds=(sp-sm)/(2*h);dd=coupling*(mp-mm)/(2*h)
                # Full dipole derivative is equivalent to the original shell
                # contraction for A, since Tr[A S]=0 for an idempotent D.
                common=2*np.einsum('ij,ij',A,dd)+np.einsum('ij,ij',pbase,ds)
                coherent=2*z*np.sqrt(omega/2)*np.einsum('ij,ij',P,dd)
                harmonic=omega*z*z/mf.mol.nelectron*np.einsum('ij,ij',P,ds)
                gold[a,k]+=common+coherent+harmonic+dexp**2/mf.mol.nelectron*np.einsum('ij,ij',D,ds)
                gnew[a,k]+=common+coherent+harmonic+dexp*np.einsum('ij,ij',P,dd)
                gfd[a,k]=(solve(plus,q,coupling).e_tot-solve(minus,q,coupling).e_tot)/(2*h)
        print('dipole=',np.einsum('ij,ji',P,mu))
        print('old forces\n',-gold)
        print('new forces\n',-gnew)
        print('finite difference\n',-gfd)
        print('old - FD forces\n',gfd-gold)
        print('max old/new errors',np.max(abs(gold-gfd)),np.max(abs(gnew-gfd)),flush=True)
        assert np.max(abs(gnew-gfd))<1e-6
        photon = -(omega**2*q+omega*coupling*np.einsum('ij,ji',P,mu))
        photon_fd = -(solve(coords,q+h,coupling).e_tot-solve(coords,q-h,coupling).e_tot)/(2*h)
        print('photon force and reference:',photon,photon_fd,'error:',abs(photon-photon_fd),flush=True)
        assert abs(photon-photon_fd)<1e-6
