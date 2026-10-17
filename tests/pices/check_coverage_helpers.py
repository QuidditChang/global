"""Exercise actual census arithmetic and the narrowly authorized metadata check."""
from pathlib import Path
import ctypes, os, random, subprocess, tempfile
ROOT=Path(__file__).resolve().parents[2]
s=(ROOT/'lib/Pices_checkpoint.c').read_text()
fn=s[s.index('static int coverage_restart_compatible('):s.index('void pices_checkpoint_check_state(')]
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp);c=p/'test.c'
 c.write_text('#include <stdlib.h>\n#include <string.h>\n#include "pices_coverage_stats.h"\n'+fn+'''
void count(long long *s,int n,int previous,int valid){pices_coverage_count(s,n,previous,valid);}
int compatible(const char *a,char *b){return coverage_restart_compatible(a,b);}
''')
 subprocess.run(['cc','-shared','-fPIC','-I'+str(ROOT/'lib'),str(c),'-o',str(p/'test.so')],check=True)
 lib=ctypes.CDLL(str(p/'test.so'));lib.count.argtypes=[ctypes.POINTER(ctypes.c_longlong),ctypes.c_int,ctypes.c_int,ctypes.c_int];lib.compatible.argtypes=[ctypes.c_char_p,ctypes.c_char_p]
 rng=random.Random(17)
 for valid in (0,1):
  values=[(rng.randrange(60),rng.randrange(60)) for _ in range(1000)];out=(ctypes.c_longlong*116)();expected=[0]*116
  for n,old in values:
   lib.count(out,n,old,valid);expected[min(n,9)]+=1
   if valid:expected[10+min(old,9)*10+min(n,9)]+=1
   for i,k in enumerate((2,4,8)):
    expected[110+2*i]+=n<k;expected[111+2*i]+=max(0,k-n)
  assert list(out)==expected
 source='bb334f35892f88b3c9a9fc5265f70b2e621049bf'
 actual=('\"solver_commit\":\"'+source+'\",\"physics\":\"abc\"').encode()
 expected=('\"solver_commit\":\"'+'1'*40+'\",\"physics\":\"abc\"').encode()
 def check(a,b):return lib.compatible(a,ctypes.create_string_buffer(b))
 os.environ.pop('PICES_RESTART_SOURCE_COMMIT',None);assert not check(actual,expected)
 os.environ['PICES_RESTART_SOURCE_COMMIT']=source;assert check(actual,expected)
 assert not check(actual.replace(b'abc',b'bad'),expected)
 assert not check(actual.replace(source.encode(),b'2'*40),expected)
 os.environ['PICES_RESTART_SOURCE_COMMIT']='2'*40;assert not check(actual,expected)
 print('PASS: census bins/transitions/deficits and strict cross-commit metadata authorization')
