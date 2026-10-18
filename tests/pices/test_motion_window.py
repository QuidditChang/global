"""Preflight format tests; synthetic files do not stand in for MPI restart tests."""
import hashlib,json,struct,tempfile,unittest
from pathlib import Path
from motion_window import prepare,parameters,SOLVER_COMMIT,standalone_runtime,standalone_mesh,restore_input_precision

class PrepareTest(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.root=Path(self.tmp.name);self.source=self.root/'source';self.source.mkdir();self.build=self.root/'build';self.build.mkdir();self.output=self.root/'output'
  exe=self.build/'PICESMotionReplay';exe.write_bytes(b'format-test-only')
  (self.build/'replay_build.json').write_text(json.dumps(dict(solver_commit=SOLVER_COMMIT,isolated_compile=True,binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())))
  (self.source/'solver_commit.txt').write_text(SOLVER_COMMIT)
  cfg='\n'.join(['nodex=2','nodey=2','nodez=2','levels=1','nproc_surf=1','nprocx=1','nprocy=1','nprocz=1','Q0=30','reference_conductivity=-1','inputdiffusivity=-1','datafile=global','datafile_old=global','datadir=DATA/%RANK','datadir_old=Restart','restart=0','solution_cycles_init=0','refstate_file=refstate_EBA.txt','coor_file=GLB.coor.global.dat'])+'\n'
  for name in ('pid1.cfg','input.cfg'):(self.source/name).write_text(cfg)
  self.template=self.root/'template.cfg';self.template.write_text(cfg)
  for name in ('refstate_EBA.txt','GLB.coor.global.dat','Inca_plateau_15Ma_v5.xy','Nazca_ridge_15Ma_v2.xy','Juan_Fernandez_Ridge_15Ma_v2.xy'):(self.source/name).write_text('fixture\n')
  d=self.source/'DATA/0';d.mkdir(parents=True)
  (d/'global.chkpt.0.0').write_bytes(struct.pack('=8i',2,2,2,1,1,1,1,0))
  (d/'global.chkpt.0.0.pices.state').write_bytes(bytes(160))
  (d/'global.chkpt.0.0.pices.manifest').write_bytes(b'x')
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
   p=parameters((self.output/('sub%d/case.cfg'%factor)).read_text());self.assertEqual(p['Q0'],'30');self.assertEqual(p['restart'],'1');self.assertEqual(p['datadir_old'],str(self.source.resolve()/'DATA/%RANK'));self.assertEqual(p['coor_file'],'GLB.coor.global.dat');self.assertNotIn('reference_conductivity',p);self.assertNotIn('inputdiffusivity',p)
  self.assertEqual(m['standalone_unset_parameters'],['inputdiffusivity','reference_conductivity'])
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

 def test_rejects_wrong_state_size_before_output(self):
  (self.source/'DATA/0/global.chkpt.0.0.pices.state').write_bytes(b'x')
  with self.assertRaisesRegex(ValueError,'state size/config mismatch'):self.run_prepare()
  self.assertFalse(self.output.exists())
 def test_rejects_wrong_mesh_before_output(self):
  p=self.source/'DATA/0/global.chkpt.0.0.pices.json';m=json.loads(p.read_text());m['local_mesh']=[3,2,2];p.write_text(json.dumps(m))
  with self.assertRaisesRegex(ValueError,'local mesh mismatch'):self.run_prepare()
  self.assertFalse(self.output.exists())

class RuntimeBridgeTest(unittest.TestCase):
 def test_explicit_values_are_preserved(self):
  for value in ('0.5859375','0','-2'):
   text='reference_conductivity='+value+'\ninputdiffusivity=-1\n'
   fixed,removed=standalone_runtime(text)
   self.assertEqual(parameters(fixed),{'reference_conductivity':value})
   self.assertEqual(removed,['inputdiffusivity'])
 def test_alias_is_preserved_when_reference_is_unset(self):
  fixed,removed=standalone_runtime('reference_conductivity=-1.0\ninputdiffusivity=0.5\nQ0=30\n')
  self.assertEqual(parameters(fixed),{'inputdiffusivity':'0.5','Q0':'30'})
  self.assertEqual(removed,['reference_conductivity'])

 def test_production_mesh_derivation(self):
  raw='nodex=129\nnodey=129\nnodez=65\nnprocx=4\nnprocy=4\nnprocz=2\nlevels=3\n'
  fixed,units,mesh=standalone_mesh(raw)
  self.assertEqual(units,{'mgunitx':8,'mgunity':8,'mgunitz':8});self.assertEqual(mesh,[33,33,33])
  with self.assertRaisesRegex(ValueError,'inconsistent mgunitx'):standalone_mesh(raw+'mgunitx=1\n')
  with self.assertRaisesRegex(ValueError,'nonintegral'):standalone_mesh(raw.replace('nodex=129','nodex=130'))
 def test_recovers_display_rounding_only(self):
  original='phase_depth=0.064354104536179568,0.1035944121801915\nrayleigh=690018487.21264744\n'
  shown='phase_depth=0.0643541,0.103594\nrayleigh=6.90018e+08\n'
  fixed,keys=restore_input_precision(shown,original)
  self.assertEqual(parameters(fixed),parameters(original));self.assertEqual(len(keys),2)
  with self.assertRaisesRegex(ValueError,'beyond %g rounding'):restore_input_precision(shown.replace('6.90018e+08','6.90019e+08'),original)

if __name__=='__main__':unittest.main()
