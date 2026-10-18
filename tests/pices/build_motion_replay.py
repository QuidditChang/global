"""Build a test-only frozen-velocity replay using the actual MPI tracer mover.
Requires objects produced by tests/cbf/build_validation.py, at the checkpoint commit.
Never advances temperature, TA, composition or Stokes during replay.
"""
import argparse,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);a=p.parse_args();b=a.build.resolve()
helper=r'''
#include "pices.h"
static void motion_replay(struct All_variables *E) {
 const char *root=getenv("PICES_REPLAY_ROOT"),*factor=getenv("PICES_REPLAY_SUBSTEPS");
 int subdivisions=factor?atoi(factor):0,step=E->monitor.solution_cycles+1,head[4],d,p;
 double clock[2];char name[2048];FILE *f;
 if(!root || (subdivisions!=1&&subdivisions!=2&&subdivisions!=4))pices_fail(E,"invalid replay controls");
 snprintf(name,sizeof(name),"%s/DATA/%d/global.motion_velocity.%d.%d",root,E->parallel.me,E->parallel.me,step);
 f=fopen(name,"rb");if(!f)pices_fail(E,"missing frozen velocity");
 if(fread(head,sizeof(int),4,f)!=4||fread(clock,sizeof(double),2,f)!=2||head[0]!=0x504d5631||head[1]!=E->parallel.me||head[2]!=step||head[3]!=E->lmesh.nno)pices_fail(E,"invalid frozen velocity header");
 if(clock[0]!=(double)E->monitor.elapsed_time)pices_fail(E,"replay checkpoint/velocity time mismatch");
 for(d=1;d<=3;d++)if(fread(E->sphere.cap[1].V[d]+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)pices_fail(E,"truncated frozen velocity");
 if(fgetc(f)!=EOF)pices_fail(E,"trailing frozen velocity data");fclose(f);
 E->monitor.solution_cycles=step;E->advection.timestep=clock[1]/subdivisions;
 for(d=0;d<subdivisions;d++)tracer_move_particles(E);
 snprintf(name,sizeof(name),"%s.replay.%d",E->control.data_file,E->parallel.me);f=fopen(name,"w");
 if(!f)pices_fail(E,"replay output open failed");
 fprintf(f,"# frozen_step=%d substeps=%d total_dt=%.17g\n",step,subdivisions,clock[1]);
 for(p=1;p<=E->trace.ntracers[1];p++)fprintf(f,"%.0f %.17g %.17g %.17g %d\n",E->trace.extraq[1][E->pices.slot+1][p],E->trace.basicq[1][3][p],E->trace.basicq[1][4][p],E->trace.basicq[1][5][p],E->trace.ielement[1][p]);
 if(fclose(f))pices_fail(E,"replay output close failed");
}
'''
s=(ROOT/'bin/Citcom.c').read_text().replace('int main(argc,argv)',helper+'\nint main(argc,argv)')
needle='      read_checkpoint(E);';assert s.count(needle)==1
s=s.replace(needle,needle+'\n      motion_replay(E);\n      MPI_Finalize();\n      return 0;')
f=b/'PICESMotionReplay.c';f.write_text(s);cc=os.environ.get('MPICC','/usr/local/bin/mpicc')
subprocess.run([cc,'-std=gnu99','-O2','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-DUSE_GZDIR','-I'+str(ROOT/'lib'),'-c',str(f),'-o',str(b/'PICESMotionReplay.o')],check=True)
sources=(ROOT/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
objects=[str(b/(Path(v).stem+'.o')) for v in sources.replace('\\','').split() if v.endswith('.c')]
subprocess.run([cc,*objects,str(b/'CitcomSFull.o'),str(b/'PICESMotionReplay.o'),'-lz','-lm','-o',str(b/'PICESMotionReplay')],check=True)
print(b/'PICESMotionReplay')
