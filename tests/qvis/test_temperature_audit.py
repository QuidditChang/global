"""Run the stage recorder on cold/negative/nonfinite nodes; verify no mutation."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]

class TemperatureAuditTest(unittest.TestCase):
    def test_absolute_kelvin_and_read_only_records(self):
        code=r'''
#include <string.h>
#include "global_defs.h"
#include "temperature_audit.h"
int main(void) {
    struct All_variables state={0}, *E=&state;
    int d; double t[4]={0,0,-0.5,NAN},td[4]={0,1,2,3};
    unsigned int flags[4]={0}; char text[8192]; size_t n;
    E->fp=tmpfile(); E->sphere.caps_per_proc=1; E->sphere.capid[1]=7;
    E->lmesh.nno=3; E->T[1]=t; E->Tdot[1]=td; E->node[1]=flags;
    E->data.Ttop=300; E->data.ref_temperature=3400;
    for(d=1;d<=3;d++) E->sx[1][d]=calloc(4,sizeof(*E->sx[1][d]));
    audit_temperature(E,"disabled",0,-1);
    if(ftell(E->fp)!=0) return 1;
    E->control.temperature_audit=1;
    audit_temperature(E,"checkpoint_loaded",0,-1);
    rewind(E->fp); n=fread(text,1,sizeof(text)-1,E->fp);text[n]=0;
    if(!strstr(text,"T_K=-1400") || !strstr(text,"negative=1 nonfinite=1")) return 2;
    if(!strstr(text,"stage=checkpoint_loaded") || !strstr(text,"cap=7 node=2")) return 3;
    if(t[1]!=0 || t[2]!=-0.5 || !isnan(t[3]) || td[2]!=2) return 4;
    return 0;
}
'''
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'test.c').write_text(code)
            subprocess.run([os.environ.get('MPICC','mpicc'),'-std=gnu99','-I'+str(ROOT/'lib'),str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
            subprocess.run([str(p/'test')],check=True)

if __name__=='__main__':unittest.main()
