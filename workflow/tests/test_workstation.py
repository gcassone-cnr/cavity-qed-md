from pathlib import Path
import os
import json
import subprocess
import sys
import tempfile
import unittest
import socket
import xml.etree.ElementTree as ET
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import et_cbohf_ipi_driver as driver
import prepare_run as prep
import launch


class WorkstationTests(unittest.TestCase):
    def test_idle_socket_eof_and_truncated_message(self):
        class FakeSocket:
            def __init__(self, chunks):self.chunks=iter(chunks)
            def recv(self, n):return next(self.chunks)
        self.assertIsNone(driver.serve(FakeSocket([b'']),None))
        with self.assertRaisesRegex(ConnectionError,'during a protocol message'):
            driver.recv_header(FakeSocket([b'S',b'']))

    def test_thermal_momenta_and_eckart_projection(self):
        names,q,p,m=prep.load_xyz(ROOT/'examples/water21.xyz',300,42)
        kinetic=np.sum(p*p/m[:,None])/2
        self.assertGreater(kinetic,0.04)
        self.assertLess(kinetic,0.16)
        projected=prep.project_external(q[:-1],p[:-1],m[:-1])
        self.assertLess(np.linalg.norm(projected.sum(axis=0)),1e-11)
        r=q[:-1]-np.average(q[:-1],axis=0,weights=m[:-1])
        self.assertLess(np.linalg.norm(np.cross(r,projected).sum(axis=0)),1e-10)

    def test_native_backend_restart_and_force_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            args=driver.build_parser().parse_args(['--et-binary',str(ROOT/'tests/harmonic_et.py'),
                 '--basis-dir',str(ROOT/'examples'),'--workdir',tmp,'--frequency','0.0187505','--coupling','0.01'])
            evaluator=driver.ETEvaluator(args,['O','H','H'])
            _,rows=driver.parse_xyz(ROOT/'examples/water.xyz')
            q=[v*driver.BOHR_PER_ANGSTROM for v in driver.flatten(rows[:-1])]+[0.5,0,0]
            _,f,md=evaluator.evaluate(q)
            self.assertFalse(md['scf_restart_loaded'])
            for i in range(10):
                plus=q.copy();minus=q.copy();plus[i]+=1e-4;minus[i]-=1e-4
                fd=-(evaluator.evaluate_energy(plus)-evaluator.evaluate_energy(minus))/2e-4
                self.assertAlmostEqual(f[i],fd,places=7)
            _,_,md=evaluator.evaluate(q)
            self.assertTrue(md['scf_restart_loaded'])
            (evaluator.restart_directory/'orbital_coefficients_qed-rhf').write_text('corrupt')
            _,_,md=evaluator.evaluate(q)
            self.assertTrue(md['scf_restart_fallback'])
            self.assertEqual((evaluator.restart_directory/'orbital_coefficients_qed-rhf').read_text(),'valid\n')

    def test_staging_rejects_eckart_fixed_dofs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);tree=ET.parse(ROOT/'examples/water_eckart.xml')
            ET.SubElement(tree.find('./system/motion'),'fixatoms_dof').text='[10,11]'
            tree.write(tmp/'bad.xml')
            with self.assertRaisesRegex(ValueError,'fixed Cartesian'):
                launch.stage(tmp/'bad.xml',tmp,'test',1)

    def test_restart_preserves_photon_and_reference(self):
        n,q,p,m=prep.load_restart(ROOT/'examples/water_eckart.xml')
        reference=q[:-1].copy()
        # Internal dilation satisfies the same frame constraints but changes q.
        com=np.average(reference,axis=0,weights=m[:-1])
        q[:-1]=com+1.02*(reference-com)
        root,_=prep.build(n,q,p,m,eckart=True,ensemble='nve',steps=5,timestep=0.1,temperature=300,stride=1,seed=1,reference=reference)
        self.assertAlmostEqual(prep.array(root.findtext('./system/beads/p'))[-3],p[-1,0])
        ref=prep.array(root.findtext('./system/motion/constrained_dynamics/constraint/values')).reshape(-1,3)
        np.testing.assert_array_equal(ref,reference)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'RESTART';ET.ElementTree(root).write(path)
            np.testing.assert_array_equal(prep.restart_reference(path,n),reference)

    def test_dft_template_is_rejected_for_cavity_md(self):
        with tempfile.TemporaryDirectory() as tmp:
            args=driver.build_parser().parse_args(['--template',str(ROOT/'examples/pbe_single_point.inp'),'--workdir',tmp])
            with self.assertRaisesRegex(ValueError,'CBO-DFT'):
                driver.ETEvaluator(args,['O','H','H'])

    def test_real_ipi_eckart_socket_and_constraints(self):
        # Real pinned i-PI integrator and sockets; harmonic test double replaces eT.
        try:
            check=socket.socket(socket.AF_UNIX,socket.SOCK_STREAM)
            check.close()
        except PermissionError:
            self.skipTest('Execution environment prohibits UNIX sockets; run this test on the workstation')
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            for name in ['water_eckart','water21_eckart','water21_nvt']:
                tree=ET.parse(ROOT/f'examples/{name}.xml');tree.find('total_steps').text='8'
                source=tmp/(name+'.xml');tree.write(source)
                env=os.environ.copy();env.update({'IPI_INPUT':str(source),'RUN_ROOT':str(tmp),
                    'RUN_TAG':name,'SOCKET_NAME':'test_'+name+'_'+str(os.getpid()),
                    'ET_BINARY':str(ROOT/'tests/harmonic_et.py'),'ET_BASIS_DIR':str(ROOT/'examples'),
                    'ET_THREADS':'1','ET_CLIENTS':'1','ET_FREQUENCY':'0.0187505','ET_COUPLING':'0.01'})
                env['ET_BACKEND']='direct'
                env.pop('ET_SCRATCH_ROOT',None)
                env.pop('ET_TEMPLATE',None)
                result=subprocess.run([sys.executable,str(ROOT/'launch.py')],env=env,capture_output=True,text=True,timeout=90)
                if result.returncode:
                    logs='\n'.join(p.read_text()[-4000:] for p in (tmp/name).glob('*.log'))
                    self.fail(result.stdout+result.stderr+logs)
                trajectory=next((tmp/name).glob('traj.positions_au*xyz'))
                lines=trajectory.read_text().splitlines();natoms=int(lines[0]);frames=[]
                for i in range(0,len(lines),natoms+2):
                    frames.append(np.array([[float(v) for v in line.split()[1:4]] for line in lines[i+2:i+2+natoms]]))
                self.assertGreaterEqual(len(frames),8)
                for frame in frames:
                    np.testing.assert_array_equal(frame[-1,1:],[0,0])
                if name.endswith('eckart'):
                    n,q,p,m=prep.load_restart(source)
                    r=q[:-1]-np.average(q[:-1],axis=0,weights=m[:-1])
                    for frame in frames:
                        delta=frame[:-1]-q[:-1]
                        self.assertLess(np.linalg.norm(np.average(delta,axis=0,weights=m[:-1])),1e-5)
                        self.assertLess(np.linalg.norm(np.cross(m[:-1,None]*r,delta).sum(axis=0)/m[:-1].sum()),1e-5)
                # Preserve the actual diagnostic evidence when requested.
                if os.environ.get('VALIDATION_LOG_DIR'):
                    out=Path(os.environ['VALIDATION_LOG_DIR']);out.mkdir(parents=True,exist_ok=True)
                    (out/(name+'_ipi.log')).write_text((tmp/name/'ipi.log').read_text())

    def test_real_ipi_internal_harmonic_eckart(self):
        import contextlib
        from ipi.engine.simulation import Simulation
        from i_pi_eckart_compat import apply_compatibility
        from ipi.utils.depend import dstrip
        apply_compatibility()
        for name in ['water_eckart','water21_eckart']:
            tree=ET.parse(ROOT/f'examples/{name}.xml');root=tree.getroot()
            root.find('total_steps').text='30'
            names,q,p,m=prep.load_restart(ROOT/f'examples/{name}.xml')
            n=q.size;H=np.eye(n)*0.001;H[-3,-3]=0.0187505**2;H[-2:,-2:]=0
            root.remove(root.find('ffsocket'))
            ff=ET.SubElement(root,'ffdebye',name='et-cbohf',pbc='false')
            ET.SubElement(ff,'hessian',shape=f'({n},{n})').text=prep.fmt(H)
            ET.SubElement(ff,'x_reference',shape=f'({n})').text=prep.fmt(q)
            out=root.find('output')
            for traj in out.findall('trajectory'):
                if traj.get('extra_type'):out.remove(traj)
            with tempfile.TemporaryDirectory() as tmp:
                previous=Path.cwd()
                try:
                    os.chdir(tmp)
                    with open('internal_ipi.log','w') as log,contextlib.redirect_stdout(log):
                        sim=Simulation.load_from_xml(ET.tostring(root,encoding='unicode'))
                        sim.run()
                        final=np.array(dstrip(sim.syslist[0].beads.q)[0]).reshape(-1,3)
                        finalp=np.array(dstrip(sim.syslist[0].beads.p)[0]).reshape(-1,3)
                        cons=sim.syslist[0].motion.constraint_list[0]
                        position_residual=float(np.linalg.norm(dstrip(cons.g)))
                        velocity_residual=float(np.linalg.norm(dstrip(cons.Dg)@(finalp[:-1]/m[:-1,None]).ravel()))
                        self.assertLess(position_residual,1e-9)
                        self.assertLess(velocity_residual,1e-10)
                        np.testing.assert_array_equal(final[-1,1:],[0,0])
                        np.testing.assert_array_equal(finalp[-1,1:],[0,0])
                    data=np.loadtxt('traj.properties')
                    self.assertLess(np.ptp(data[:,2]),1e-4) # eV, short harmonic NVE test only
                    if os.environ.get('VALIDATION_LOG_DIR'):
                        outpath=Path(os.environ['VALIDATION_LOG_DIR']);outpath.mkdir(parents=True,exist_ok=True)
                        (outpath/(name+'_internal_ipi.log')).write_text(Path('internal_ipi.log').read_text())
                        (outpath/(name+'_harmonic_metrics.json')).write_text(json.dumps({
                            'model':'i-PI internal harmonic force, not eT', 'steps':30,
                            'timestep_fs':0.1,'eckart_position_residual_norm':position_residual,
                            'eckart_velocity_residual_norm':velocity_residual,
                            'energy_peak_to_peak_eV':float(np.ptp(data[:,2]))},indent=2)+'\n')
                finally:
                    os.chdir(previous)


if __name__=='__main__':
    unittest.main()
