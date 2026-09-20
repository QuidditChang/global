"""Check native CBF equations, face integrals and expected rank files.
Usage: python verify_native_outputs.py RUN_DIR [--boundary-ranks 12] [--step 1]
"""
import argparse
from pathlib import Path
import re
import math

def verify(root, step, ranks):
    for prefix, sign in [('botm', -1), ('surf', 1)]:
        files = sorted(root.rglob('q.%s.*.%d' % (prefix,step)))
        assert len(files)==ranks, (prefix,len(files),ranks)
        total=0.; nodes_count=0; faces_count=0; expected=None
        for path in files:
            nodes={};faces=[];meta={}
            for line in path.read_text().splitlines():
                if line.startswith('#'):
                    meta.update(re.findall(r'(\w+)=([^\s]+)',line))
                elif line.startswith('N '):
                    parts=line.split(); key=tuple(map(int,parts[1:3])); values=list(map(float,parts[3:]))
                    assert all(map(math.isfinite,values)),path
                    nodes[key]=values
                elif line.startswith('F '): faces.append(line.split())
            length=float(meta['length_scale_m']);scale=float(meta['k0_W_m_K'])*float(meta['deltaT_K'])/length
            target=float(meta['global_heat_W'])
            if expected is None: expected=target
            assert math.isclose(expected,target,rel_tol=1e-14,abs_tol=1e-8)
            for values in nodes.values():
                q,b,mass=values[-3:];assert mass>0
                assert math.isclose(q,sign*scale*b/mass,rel_tol=1e-13,abs_tol=1e-13)
            for f in faces:
                cap=int(f[1]); ids=list(map(int,f[3:7]));weights=list(map(float,f[7:11]))
                assert all(w>0 for w in weights)
                total += sum(nodes[(cap,n)][-3]*w*length**2 for n,w in zip(ids,weights))
            nodes_count+=len(nodes);faces_count+=len(faces)
        assert math.isclose(total,expected,rel_tol=1e-12,abs_tol=1.),(total,expected)
        print('%s step=%d ranks=%d nodes=%d faces=%d power_W=%.12g PASS' %
              (prefix,step,ranks,nodes_count,faces_count,total))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path)
    p.add_argument('--boundary-ranks',type=int,default=12);p.add_argument('--step',type=int,default=1)
    a=p.parse_args();verify(a.root,a.step,a.boundary_ranks)
