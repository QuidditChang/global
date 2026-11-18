"""Compile production 7/71/72 kernels and verify physical contracts.

Usage: python3 tests/test_rheol7_nuref.py WORKSPACE_OUTPUT
The independently fitted 65-node M2-A archive is mode 72's compatibility oracle.
Dividing the saved viscosity by its saved correction factor gives mode 7's oracle.
"""
from pathlib import Path
import math
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
out = Path(sys.argv[1]).resolve()
source = r'''
#include <assert.h>
#include <stdio.h>
#include "Steinberger_nuref.h"
static int close_to(double a,double b) {return fabs(a/b-1.)<2.e-13;}
static void contracts(void) {
 double q[4],a[4],b[4],c[4];
 const int modes[3]={7,71,72};
 const double depths[]={0.,200.,410.-1.e-7,410.,520.-1.e-7,520.,660.-1.e-7,660.,1000.,2891.};
 int i,j;
 steinberger_defaults(q,a,b,c);
 for(i=0;i<3;i++) {
  assert(steinberger_is_family(modes[i]));
  assert(steinberger_validate_correction(modes[i],q,a,b,c));
 }
 assert(!steinberger_is_family(70));assert(!steinberger_is_family(73));
 assert(!steinberger_validate_correction(73,q,a,b,c));
 /* Identity correction in both parameterizations. */
 for(i=0;i<4;i++) {q[i]=a[i]=1.;b[i]=c[i]=0.;}
 for(j=0;j<10;j++) {
  double z=depths[j],base=steinberger_nuref(z,2000.,1.e21);
  for(i=0;i<3;i++) {
   assert(close_to(steinberger_nuref_mode(modes[i],z,2000.,1.e21,q,a,b,c),base));
   assert(close_to(steinberger_viscosity_mode(modes[i],z,2000.,1800.,1.e21,.5,q,a,b,c),
                   steinberger_viscosity(z,2000.,1800.,1.e21,.5)));
  }
 }
 /* Custom corrections preserve common cold/hot response; deeper layer owns boundaries. */
 for(i=0;i<4;i++) {q[i]=2.+i;a[i]=1.+i;b[i]=.5;c[i]=.25;}
 for(j=0;j<10;j++) {
  double z=depths[j],bounds[5]={0.,410.,520.,660.,2891.};
  int layer=z<410.?0:(z<520.?1:(z<660.?2:3));
  double u=(z-bounds[layer])/(bounds[layer+1]-bounds[layer]);
  double expected=1.+layer+.5*u+.25*u*u;
  assert(close_to(steinberger_factor(71,z,q,a,b,c),2.+layer));
  assert(close_to(steinberger_factor(72,z,q,a,b,c),expected));
  for(i=0;i<3;i++) {
   double f=i==0?1.:(i==1?2.+layer:expected);
   double ref=steinberger_nuref_mode(modes[i],z,2000.,1.e21,q,a,b,c);
   assert(close_to(ref,steinberger_nuref(z,2000.,1.e21)*f));
   assert(close_to(2.*steinberger_nuref_mode(modes[i],z,2000.,2.e21,q,a,b,c),ref));
   assert(close_to(steinberger_viscosity_mode(modes[i],z,2000.,1800.,1.e21,.5,q,a,b,c),
                   steinberger_viscosity(z,2000.,1800.,1.e21,.5)*f));
   assert(close_to(steinberger_viscosity_mode(modes[i],z,2000.,2200.,1.e21,.5,q,a,b,c),
                   steinberger_viscosity(z,2000.,2200.,1.e21,.5)*f));
  }
 }
 for(i=0;i<3;i++) {
  double z=1000.,r=2100.,e=1.e21;
  double ref=steinberger_nuref_mode(modes[i],z,r,e,q,a,b,c);
  double full=steinberger_viscosity_mode(modes[i],z,r,1900.,e,1.,q,a,b,c);
  double half=steinberger_viscosity_mode(modes[i],z,r,1900.,e,.5,q,a,b,c);
  double quarter=steinberger_viscosity_mode(modes[i],z,r,1900.,e,.25,q,a,b,c);
  assert(steinberger_viscosity_mode(modes[i],z,r,1900.,e,0.,q,a,b,c)==ref);
  assert(fabs(log(half/ref)*2.-log(full/ref))<1.e-12);
  assert(fabs(log(quarter/ref)*4.-log(full/ref))<1.e-12);
  assert(steinberger_viscosity_mode(modes[i],z,r,2300.,e,0.,q,a,b,c)==steinberger_viscosity_mode(modes[i],z,r,2300.,e,1.,q,a,b,c));
  assert(steinberger_viscosity_mode(modes[i],z,r,1900.,e,-1.,q,a,b,c)<0.);
  assert(steinberger_viscosity_mode(modes[i],z,r,1900.,e,NAN,q,a,b,c)<0.);
  assert(steinberger_viscosity_mode(modes[i],z,r,NAN,e,1.,q,a,b,c)<0.);
  assert(steinberger_nuref_mode(modes[i],-1.,r,e,q,a,b,c)<0.);
  assert(steinberger_nuref_mode(modes[i],2892.,r,e,q,a,b,c)<0.);
  assert(steinberger_nuref_mode(modes[i],NAN,r,e,q,a,b,c)<0.);
  assert(steinberger_nuref_mode(modes[i],z,0.,e,q,a,b,c)<0.);
  assert(steinberger_nuref_mode(modes[i],z,r,0.,q,a,b,c)<0.);
 }
 /* Whole-interval positivity: endpoints alone are insufficient. */
 steinberger_defaults(q,a,b,c);
 q[2]=0.;assert(!steinberger_validate_correction(71,q,a,b,c));
 q[2]=-1.;assert(!steinberger_validate_correction(71,q,a,b,c));
 q[2]=NAN;assert(!steinberger_validate_correction(71,q,a,b,c));
 q[2]=1.;a[0]=1.;b[0]=-8.;c[0]=8.;
 assert(!steinberger_validate_correction(72,q,a,b,c));
 b[0]=-4.;c[0]=4.;assert(!steinberger_validate_correction(72,q,a,b,c));
 a[0]=1.01;assert(steinberger_validate_correction(72,q,a,b,c));
 b[0]=INFINITY;assert(!steinberger_validate_correction(72,q,a,b,c));
 b[0]=0.;c[0]=NAN;assert(!steinberger_validate_correction(72,q,a,b,c));
 c[0]=-2.;assert(!steinberger_validate_correction(72,q,a,b,c));
}
int main(void) {
 double q[4],a[4],b[4],c[4],z,t,e;int mode;
 contracts();steinberger_defaults(q,a,b,c);
 while(scanf("%d %lf %lf %lf",&mode,&z,&t,&e)==4)
  printf("%.17g\n",steinberger_nuref_mode(mode,z,t,e,q,a,b,c));
 return 0;
}
'''
nodes = [list(map(float,line.split())) for line in
         (out/'steinberger_M2_polynomial_factors/M2_A_65_nodes.txt').read_text().splitlines()
         if line and not line.startswith('#')]
with tempfile.TemporaryDirectory() as temporary:
    path = Path(temporary)
    (path/'test.c').write_text(source)
    subprocess.run(['cc','-std=c99','-Wall','-Wextra','-Werror','-I'+str(root/'lib'),
                    str(path/'test.c'),'-lm','-o',str(path/'test')],check=True)
    def evaluate(mode):
        data=''.join(f'{mode} {r[2]:.17g} {r[4]:.17g} 1e21\n' for r in nodes)
        return list(map(float,subprocess.check_output([str(path/'test')],input=data,text=True).split()))
    actual=evaluate(72)
    error=max(abs(a/r[7]-1) for a,r in zip(actual,nodes));assert error<1e-9,error
    baseline=evaluate(7)
    base_error=max(abs(a/(r[7]/r[5])-1) for a,r in zip(baseline,nodes));assert base_error<1e-9,base_error
    assert all(math.isfinite(v) and v>0 for v in actual+baseline)
    assert all(abs(a/b-1)<2e-13 for a,b in zip(evaluate(71),baseline))
print(f'PASS: mode72 archived65 maximum relative error {error:.3e}; '
      f'mode7 baseline oracle error {base_error:.3e}; mode71 identity.')
print('PASS: all-mode thermal response/refvisc scaling; custom factors; '
      'deeper-side boundaries; whole-interval positivity; invalid inputs.')
