"""Production checkpoint clock reader: precision, early init and invalid headers."""
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from test_cbf_kernel import function, ROOT

class RestartClock(unittest.TestCase):
    def test_binary_header_and_initialization_cycle(self):
        code=function((ROOT/'lib/Checkpoints.c').read_text(),'void read_checkpoint_initial_time(')
        prefix=r'''
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_MIN 1
#define MPI_MAX 2
static void MPI_Allreduce(void *a,void *b,int n,int t,int op,int world) {
 memcpy(b,a,n*(t==MPI_FLOAT ? sizeof(float):sizeof(int)));
}
struct All_variables {
 struct {char old_P_file[200];float start_age;} control;
 struct {int solution_cycles_init,solution_cycles;float elapsed_time;} monitor;
 struct {int me,nprocx,nprocy,nprocz,world;} parallel;
 struct {int nox,noy,noz;} lmesh;
 struct {int caps_per_proc;} sphere;
 FILE *fp;
};
void parallel_process_termination(void) {exit(42);}
'''
        main=r'''
int main(int argc,char **argv) {
 struct All_variables E={0};
 snprintf(E.control.old_P_file,200,"%s/global",argv[1]);
 E.control.start_age=123; E.monitor.solution_cycles_init=13600;
 E.lmesh.nox=E.lmesh.noy=E.lmesh.noz=33;
 E.parallel.nprocx=E.parallel.nprocy=4;E.parallel.nprocz=2;
 E.sphere.caps_per_proc=1;E.fp=stdout;
 read_checkpoint_initial_time(&E);
 if(E.monitor.solution_cycles!=0)return 3;
 if(E.monitor.elapsed_time!=(float)0.000193788626348)return 4;
 if(E.control.start_age!=(float)249.9)return 5;
 return 0;
}
'''
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);src=root/'clock.c';src.write_text(prefix+code+main)
            exe=root/'clock';subprocess.run(['cc','-std=gnu99',str(src),'-o',str(exe)],check=True)
            path=root/'global.chkpt.0.13600'
            valid=struct.pack('=8i3f',33,33,33,4,4,2,1,13600,.000193788626348,1.e-8,249.9)
            path.write_bytes(valid)
            subprocess.run([str(exe),d],check=True,capture_output=True)
            for bad in (valid[:20],struct.pack('=8i3f',33,33,33,4,4,2,1,7,.1,1.e-8,249.9),
                        struct.pack('=8i3f',33,33,33,4,4,2,1,13600,float('nan'),1.e-8,249.9)):
                path.write_bytes(bad)
                result=subprocess.run([str(exe),d],capture_output=True)
                self.assertEqual(result.returncode,42)

if __name__=='__main__':unittest.main()
