"""Check the production C sampler against independent area quadrature.

The asymmetric triangles expose incorrect plane-area mixture weights and a
missing gnomonic Jacobian. Radial moments expose uniform-r instead of r^3.
"""
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HARNESS = r'''
#include <stdio.h>
#include "tracer_seed_math.h"
int main(void) {
 const double uv[4][2]={{-.8,-.6},{.7,-.4},{.9,.8},{-.2,.5}};
 double area[2]={tracer_seed_area(uv[0],uv[1],uv[2]),tracer_seed_area(uv[0],uv[2],uv[3])};
 uint64_t state=20261114;long n=500000,i;double su=0,sv=0,sq=0,sq2=0,suq=0;
 for(i=0;i<n;i++) {
  double u,v,r,q;
  if(!tracer_seed_direction(&state,uv,area,&u,&v))return 1;
  r=tracer_seed_radius(&state,.55,1);q=(r*r*r-.55*.55*.55)/(1-.55*.55*.55);
  if(!(q>0&&q<1))return 2;
  su+=u;sv+=v;sq+=q;sq2+=q*q;suq+=u*q;
 }
 printf("{\"u\":%.17g,\"v\":%.17g,\"q\":%.17g,\"q2\":%.17g,\"uq\":%.17g}\n",su/n,sv/n,sq/n,sq2/n,suq/n);
 return 0;
}
'''


def main():
    with tempfile.TemporaryDirectory(prefix='particle-seed-sampling-') as name:
        temp = Path(name)
        (temp / 'check.c').write_text(HARNESS)
        subprocess.run([os.environ.get('CC', 'cc'), '-std=c99', '-O2',
                        '-I' + str(ROOT / 'lib'), str(temp / 'check.c'),
                        '-lm', '-o', str(temp / 'check')], check=True)
        got = json.loads(subprocess.check_output([str(temp / 'check')], text=True))
    vertices = ((-.8, -.6), (.7, -.4), (.9, .8), (-.2, .5))
    mass = u_moment = v_moment = 0.0
    resolution = 600
    for index in ((0, 1, 2), (0, 2, 3)):
        a, b, c = (vertices[i] for i in index)
        area = abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])) / 2
        for i in range(resolution):
            s = math.sqrt((i + .5) / resolution)
            for j in range(resolution):
                t = (j + .5) / resolution
                u = (1-s)*a[0] + s*(1-t)*b[0] + s*t*c[0]
                v = (1-s)*a[1] + s*(1-t)*b[1] + s*t*c[1]
                weight = area / (1 + u*u + v*v)**1.5
                mass += weight
                u_moment += weight * u
                v_moment += weight * v
    expected = dict(u=u_moment / mass, v=v_moment / mass, q=.5, q2=1/3)
    for key, value in expected.items():
        assert abs(got[key] - value) < .003, (key, got[key], value)
    assert abs(got['uq'] - got['u'] * got['q']) < .0015, got
    print(json.dumps(dict(status='PASS', sampled=got, expected=expected,
                         samples=500000, checks=['spherical-area Jacobian and triangle mixture',
                                                 'radial volume moments', 'radial/angular independence']), indent=2))


if __name__ == '__main__':
    main()
