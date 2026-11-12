#!/usr/bin/env python3
"""Compile the production analytic map and check its conservation identities."""
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class VolumeControlMath(unittest.TestCase):
    def test_compact_vortex_and_exact_volume_preservation(self):
        compiler = shutil.which("cc")
        self.assertIsNotNone(compiler, "a C compiler is required")
        source = r'''
#include <stdio.h>
#include "pices_volume_control_math.h"
static double norm(const double x[3]) {
    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}
static double determinant(double a[3][3]) {
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
int main(void) {
    struct Pices_control_flow f={{.64/sqrt(2.04),.64/sqrt(2.04),.128/sqrt(2.04)},.08,acos(-1.0)/4e-6};
    double t=4e-6,maxdet=0,maxdiv=0,maxroundtrip=0,h=1e-7;
    int i,j,k;
    for(i=0;i<1600;i++) {
        double x[3],y[3],back[3],jac[3][3],div=0,er=0,rad0=0,rad1=0;
        /* Irrational quasi-uniform points include both sides of the support. */
        for(j=0;j<3;j++)x[j]=f.center[j]+.2*(fmod((i+1)*sqrt(j+2.5),1.0)-.5);
        pices_control_exact(&f,x,t,y);pices_control_exact(&f,y,-t,back);
        for(j=0;j<3;j++) {
            double xp[3],xm[3],yp[3],ym[3],up[3],um[3];
            for(k=0;k<3;k++){xp[k]=xm[k]=x[k];}
            xp[j]+=h;xm[j]-=h;
            pices_control_exact(&f,xp,t,yp);pices_control_exact(&f,xm,t,ym);
            pices_control_velocity(&f,xp,up);pices_control_velocity(&f,xm,um);
            div+=(up[j]-um[j])/(2*h);
            for(k=0;k<3;k++)jac[k][j]=(yp[k]-ym[k])/(2*h);
            er=fmax(er,fabs(back[j]-x[j]));
            rad0+=(x[j]-f.center[j])*(x[j]-f.center[j]);
            rad1+=(y[j]-f.center[j])*(y[j]-f.center[j]);
        }
        maxdet=fmax(maxdet,fabs(determinant(jac)-1));
        maxdiv=fmax(maxdiv,fabs(div/f.omega));maxroundtrip=fmax(maxroundtrip,er);
        if(fabs(rad0-rad1)>1e-15)return 1;
    }
    if(maxdet>2e-8 || maxdiv>2e-8 || maxroundtrip>3e-15)return 2;
    {
        double x[3],u[3],y[3],bound=16*f.omega*f.radius/(25*sqrt(5.0));
        for(j=0;j<3;j++)x[j]=f.center[j];x[0]+=f.radius/sqrt(5.0);
        pices_control_velocity(&f,x,u);
        if(fabs(norm(u)/bound-1)>1e-12)return 3;
        x[0]=f.center[0]+1.01*f.radius;pices_control_exact(&f,x,t,y);
        pices_control_velocity(&f,x,u);
        for(j=0;j<3;j++)if(x[j]!=y[j] || u[j]!=0)return 4;
        for(j=0;j<3;j++)x[j]=f.center[j]*.575/.64;
        pices_control_exact(&f,x,t,y);
        if(norm(y)-norm(x)<.001)return 5;
    }
    printf("max |det J-1|=%.3g, max |div u|/omega=%.3g, inverse error=%.3g\n",maxdet,maxdiv,maxroundtrip);
    return 0;
}
'''
        with tempfile.TemporaryDirectory(prefix="pices-control-math-") as td:
            path = pathlib.Path(td)
            (path / "test.c").write_text(source)
            subprocess.run([compiler, "-std=c99", "-O2", "-I", str(ROOT / "lib"),
                            str(path / "test.c"), "-lm", "-o", str(path / "test")], check=True)
            result = subprocess.run([str(path / "test")], check=True, text=True, capture_output=True)
            print(result.stdout, end="")


if __name__ == "__main__":
    unittest.main()
