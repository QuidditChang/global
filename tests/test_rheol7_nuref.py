"""Compile the production kernel and compare with the independently fitted profile.
Run: python3 tests/test_rheol7_nuref.py /path/to/CitcomS/output
"""
import pathlib, subprocess, tempfile, sys, math
root=pathlib.Path(__file__).resolve().parents[1]
out=pathlib.Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as tmp:
 p=pathlib.Path(tmp)
 (p/'test.c').write_text('#include <stdio.h>\n#include "Steinberger_nuref.h"\nint main(void) {double z,t,e; while(scanf("%lf %lf %lf",&z,&t,&e)==3) printf("%.17g\\n",steinberger_viscosity(z,t,t,e,0.5)); return 0;}\n')
 subprocess.run(['cc','-std=c99','-Wall','-Wextra','-Werror','-I'+str(root/'lib'),str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
 def run(rows):
  data=''.join(' '.join(map(str,r))+'\n' for r in rows)
  return list(map(float,subprocess.check_output([str(p/'test')],input=data,text=True).split()))
 nodes=[list(map(float,l.split())) for l in (out/'steinberger_M2_polynomial_factors/M2_A_65_nodes.txt').read_text().splitlines() if l and not l.startswith('#')]
 actual=run([(r[2],r[4],1e21) for r in nodes])
 error=max(abs(a/r[7]-1) for a,r in zip(actual,nodes));assert error<1e-9,error
 scaled=run([(r[2],r[4],2e21) for r in nodes]);assert all(abs(2*b/a-1)<1e-13 for a,b in zip(actual,scaled))
 # Equal temperatures isolate the intended segment assignment from Tref jumps.
 boundaries=run([(z,2000,1e21) for z in [410-1e-7,410,520-1e-7,520,660-1e-7,660]])
 assert all(math.isfinite(a) and a>0 for a in boundaries)
 assert boundaries[0]!=boundaries[1] and boundaries[2]!=boundaries[3] and boundaries[4]!=boundaries[5]
 bad=run([(-1,2000,1e21),(2892,2000,1e21),(0,0,1e21),(0,2000,0),(float('nan'),2000,1e21)])
 assert bad==[-1]*5
 print(f'PASS: 65-node maximum relative error {error:.3e}; refvisc scaling; phase jumps; invalid-input rejection.')

# Compile parameter checks against the actual production function.
with tempfile.TemporaryDirectory() as tmp:
 p=pathlib.Path(tmp)
 (p/'cold_scale.c').write_text(r'''#include <assert.h>
#include "Steinberger_nuref.h"
int main(void) {
 double z=1000.,r=2100.,e=1e21;
 double ref=steinberger_nuref(z,r,e);
 double half=steinberger_viscosity(z,r,1900.,e,.5);
 double quarter=steinberger_viscosity(z,r,1900.,e,.25);
 double full=steinberger_viscosity(z,r,1900.,e,1.);
 assert(steinberger_viscosity(z,r,1900.,e,0.)==ref);
 assert(fabs(log(quarter/ref)*4-log(full/ref))<1e-12);
 assert(fabs(log(half/ref)*2-log(full/ref))<1e-12);
 assert(steinberger_viscosity(z,r,2300.,e,0.)==steinberger_viscosity(z,r,2300.,e,1.));
 assert(steinberger_viscosity(z,r,1900.,e,-1.)<0.);
 assert(steinberger_viscosity(z,r,1900.,e,NAN)<0.);
 return 0;
}
''')
 subprocess.run(['cc','-std=c99','-Wall','-Wextra','-Werror','-I'+str(root/'lib'),str(p/'cold_scale.c'),'-lm','-o',str(p/'cold_scale')],check=True)
 subprocess.run([str(p/'cold_scale')],check=True)
 print('PASS: cold_scale=0/0.25/0.5/1 cold response, hot-side independence, invalid cold_scale rejection.')
