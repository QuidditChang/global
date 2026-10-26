"""Compile the exact old-checkpoint authorization comparator; no hash bypass."""
import ctypes
import os
from pathlib import Path
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parents[2]
s=(ROOT/'lib/Pices_checkpoint.c').read_text()
fn=s[s.index('static int uncovered_restart_compatible('):s.index('void pices_checkpoint_check_state(')]
with tempfile.TemporaryDirectory() as folder:
 p=Path(folder)
 (p/'test.c').write_text('#include <stdlib.h>\n#include <string.h>\n'+fn+'\nint check(const char *a,char *b){return uncovered_restart_compatible(a,b);}\n')
 subprocess.run(['cc','-shared','-fPIC',str(p/'test.c'),'-o',str(p/'test.so')],check=True)
 lib=ctypes.CDLL(str(p/'test.so'));lib.check.argtypes=[ctypes.c_char_p,ctypes.c_char_p]
 def check(a,b):return lib.check(a.encode(),ctypes.create_string_buffer(b.encode()))
 def meta(commit):return '{"solver_commit":"'+commit+'","physics_mesh_sha256":"PHYSICS","checkpoint_sha256":"PAYLOAD","accepted_velocity_sha256":"STATE","time":1.25}'
 current=meta('f'*40)
 for source in ['43d530bde585b49e3dec2ab22cdb2d62e89f9af1','6a8401420be67b7f280b885e185bd049e711c4c7']:
  actual=meta(source)
  os.environ.pop('PICES_RESTART_SOURCE_COMMIT',None);assert not check(actual,current)
  os.environ['PICES_RESTART_SOURCE_COMMIT']=source;assert check(actual,current)
  for word in ['PHYSICS','PAYLOAD','STATE','1.25']:
   assert not check(actual.replace(word,'BAD'),current),word
  assert not check(meta('a'*40),current)
 os.environ['PICES_RESTART_SOURCE_COMMIT']='a'*40;assert not check(meta('a'*40),current)
 # New files encode the new policy; only the explicit legacy verification path
 # recomputes the old policy fingerprint, never substitutes a stored hash.
 assert 'if(!legacy_coverage && !E->pices.consistent_projection)' in s
 assert 'metadata_policy(E,payload,state,out,0)' in s
 assert 'metadata_policy(E,path,state,expected,1)' in s
 print('PASS: explicit audited transition; mismatched physics, payload, velocity, time and unaudited commits rejected')
