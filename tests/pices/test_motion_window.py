"""Preflight format tests; synthetic files do not stand in for MPI restart tests."""
import hashlib,json,struct,tempfile,unittest
from pathlib import Path
from motion_window import prepare,parameters,SOLVER_COMMIT

class PrepareTest(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.root=Path(self.tmp.name);self.source=self.root/'source';self.source.mkdir();self.build=self.root/'build';self.build.mkdir();self.output=self.root/'output'
  exe=self.build/'PICESMotionReplay';exe.write_bytes(b'format-test-only')
  (self.build/'replay_build.json').write_text(json.dumps(dict(solver_commit=SOLVER_COMMIT,isolated_compile=True,binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())))
  (self.source/'solver_commit.txt').write_text(SOLVER_COMMIT)
  cfg='\n'.join(['nproc_surf=1','nprocx=1','nprocy=1','nprocz=1','Q0=30','datafile=global','datafile_old=global','datadir=DATA/%RANK','datadir_old=Restart','restart=0','solution_cycles_init=0','refstate_file=refstate_EBA.txt','coor_file=GLB.coor.global.dat'])+'\n'
  for name in ('pid1.cfg','input.cfg'):(self.source/name).write_text(cfg)
  self.template=self.root/'template.cfg';self.template.write_text(cfg)
  for name in ('refstate_EBA.txt','GLB.coor.global.dat','Inca_plateau_15Ma_v5.xy','Nazca_ridge_15Ma_v2.xy','Juan_Fernandez_Ridge_15Ma_v2.xy'):(self.source/name).write_text('fixture\n')
  d=self.source/'DATA/0';d.mkdir(parents=True)
  for suffix in ('','.pices.state','.pices.manifest'):(d/('global.chkpt.0.0'+suffix)).write_bytes(b'x')
  (d/'global.chkpt.0.0.pices.json').write_text(json.dumps(dict(solver_commit=SOLVER_COMMIT,step=0,time=0.,local_mesh=[2,2,2],extraq=[{'name':'motion_id'},{'name':'motion_selected'}])))
  clock=0.;dt=struct.unpack('=f',struct.pack('=f',1e-8))[0]
  for step in (1,2):
   (d/('global.motion_velocity.0.%d'%step)).write_bytes(struct.pack('=4i2d24f',0x504d5631,0,step,8,clock,dt,*([0.]*24)))
   clock=struct.unpack('=f',struct.pack('=f',clock+dt))[0]
  (d/'global.motion.0').write_text(''.join('PATH %d 2 1 1 1e-8 .7 0 0 .7 0 0 .7 0 0 0 0 0\nHOST %d 2 1 1 1\nCELL %d 1 1\n'%(s,s,s+1) for s in (1,2)))
  (d/'global.log').write_text('PICES_MOTION step=1 stage=1 count=1 hash_sum=1 hash_xor=1\n'+''.join('PICES_COVERAGE_STATS step=%d histogram=0,1,0,0,0,0,0,0,0,0\n'%s for s in (1,2)))
 def run_prepare(self):prepare(self.source,self.output,self.template,self.build,first=1,last=2,ranks=1)
 def test_reuses_normalized_parameters_and_tracks_source(self):
  self.run_prepare();m=json.loads((self.output/'replay_manifest.json').read_text());self.assertEqual(m['tracked'],1);self.assertEqual(m['low_cells']['2'],[[0,1,1]])
  for factor in (1,2,4):
   p=parameters((self.output/('sub%d/case.cfg'%factor)).read_text());self.assertEqual(p['Q0'],'30');self.assertEqual(p['restart'],'1');self.assertEqual(p['datadir_old'],str(self.source.resolve()/'DATA/%RANK'));self.assertEqual(p['coor_file'],'GLB.coor.global.dat')
 def test_rejects_changed_physics_before_creating_output(self):
  self.template.write_text(self.template.read_text().replace('Q0=30','Q0=31'))
  with self.assertRaisesRegex(ValueError,'changed physics'):self.run_prepare()
  self.assertFalse(self.output.exists())
 def test_rejects_discontinuous_clock(self):
  p=self.source/'DATA/0/global.motion_velocity.0.2';raw=bytearray(p.read_bytes());struct.pack_into('=d',raw,16,2.);p.write_bytes(raw)
  with self.assertRaisesRegex(ValueError,'discontinuous'):self.run_prepare()
 def test_rejects_changed_binary(self):
  (self.build/'PICESMotionReplay').write_bytes(b'changed')
  with self.assertRaisesRegex(ValueError,'binary/build manifest'):self.run_prepare()

if __name__=='__main__':unittest.main()
