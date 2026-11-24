"""Test strict source-only compatibility and the pinned restart preparation."""
from pathlib import Path
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
# Local test points at the connected scripts checkout; production uses /work.
SCRIPTS=ROOT.parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('resume_host_fix',ROOT/'tools/1p8ga_resume/prepare.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
CFG='''[CitcomS]
steps=300000
[CitcomS.job]
walltime=720 *hour
[CitcomS.controller]
checkpointFrequency=50
[CitcomS.solver]
datafile=global
datadir=old
 datadir_old=old
[CitcomS.solver.ic]
restart=on
solution_cycles_init=7700
zero_elapsed_time=off
[CitcomS.solver.tsolver]
energy_solver=pices
pices_cfl=0.5
[CitcomS.solver.mesher]
nproc_surf=1
nprocx=1
nprocy=1
nprocz=1
'''.replace(' datadir_old','datadir_old')
class Restart(unittest.TestCase):
 def test_exact_step_and_original_physics(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)/'run';root.mkdir();storage=Path(tmp)/'storage';storage.mkdir()
   for d in ['DATA/0','PostProc/Online/ready','Runtime_fixed/install']: (root/d).mkdir(parents=True)
   (root/'input.cfg').write_text(CFG)
   (root/'DATA/0/global.coord.0').write_text('coordinates')
   (root/'DATA/0/global.chkpt.0.39400.pices.json').write_text(json.dumps(dict(solver_commit=m.SOURCE,magic='CITCOMS_EBA_PICES')))
   (root/'PostProc/Online/state.json').write_text(json.dumps(dict(done=[],failed={})))
   exe=root/'Runtime_fixed/install/exe';exe.write_text('fixed')
   (root/'Runtime_fixed/manifest.json').write_text(json.dumps(dict(solver_commit='a'*40,copies=[dict(destination=str(exe.parent),sha256={'exe':m.digest(exe)})])))
   with patch.object(m.Cleanup,'validate_checkpoint') as check:
    r=m.prepare(root,root/'input.cfg',storage,'999',root/'Runtime_fixed',m.SOURCE)
   check.assert_called_once_with(39400)
   self.assertEqual(r['step'],39400)
   text=(root/'input.cfg').read_text()
   self.assertIn('solution_cycles_init = 39400',text)
   self.assertIn('pices_cfl=0.5',text)
   self.assertIn('steps = 300000',text)
   self.assertEqual((root/'PostProc/Online/resume-999/input.cfg').read_text(),CFG)
 def test_c_checkpoint_identity_gate(self):
  src=(ROOT/'lib/Pices_checkpoint.c').read_text()
  fn=src[src.index('static int coverage_restart_compatible('):src.index('/* Explicit, audited transition')]
  c='''#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
'''+fn+r'''
int main(void){
 const char *source="0608be89a18479a0d22374eb9b395be940c49d86";
 char a[256],b[256],base[256];
 sprintf(a,"{\"solver_commit\":\"%s\",\"physics\":\"same\"}",source);
 sprintf(base,"{\"solver_commit\":\"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\",\"physics\":\"same\"}");
 unsetenv("PICES_RESTART_SOURCE_COMMIT");strcpy(b,base);assert(!coverage_restart_compatible(a,b));
 setenv("PICES_RESTART_SOURCE_COMMIT",source,1);strcpy(b,base);assert(coverage_restart_compatible(a,b));
 strcpy(b,base);b[strlen(b)-4]='X';assert(!coverage_restart_compatible(a,b));
 setenv("PICES_RESTART_SOURCE_COMMIT","bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",1);
 strcpy(b,base);assert(!coverage_restart_compatible(a,b));
 return 0;
}
'''
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'gate.c').write_text(c)
   subprocess.run(['cc',str(p/'gate.c'),'-o',str(p/'gate')],check=True)
   subprocess.run([str(p/'gate')],check=True)
if __name__=='__main__':unittest.main()
