/* Opt-in, observational particle-motion audit. No state feedback to physics.
 * Two extraq slots travel through the existing MPI/checkpoint machinery.
 * ID = initial rank * 2^32 + initial index (exact in a double for this layout).
 * Native binary replay dumps carry a format marker; do not add to small bundles.
 */
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "global_defs.h"
#include "pices.h"
static FILE *paths;
static int *counts,*selected;
static double *mask;
static unsigned long long reference[3],counters[3];
static double maxima[2];
int pices_motion_enabled(struct All_variables *E) {
 static int enabled=-1;
 if(enabled<0){const char *v=getenv("PICES_MOTION_AUDIT");enabled=v&&!strcmp(v,"1");}
 return E->pices.enabled && enabled;
}
static unsigned long long mix(unsigned long long x) {
 x^=x>>30;x*=UINT64_C(0xbf58476d1ce4e5b9);x^=x>>27;
 x*=UINT64_C(0x94d049bb133111eb);return x^(x>>31);
}
static void signatures(struct All_variables *E,unsigned long long out[3]) {
 unsigned long long local[3]={0,0,0};int p;
 for(p=1;p<=E->trace.ntracers[1];p++) {
  double v=E->trace.extraq[1][E->pices.slot+1][p];
  if(!isfinite(v)||v<1||v>=9007199254740992.0||v!=floor(v))pices_fail(E,"invalid motion ID");
  unsigned long long h=mix((unsigned long long)v);local[0]++;local[1]+=h;local[2]^=h;
 }
 MPI_Allreduce(local,out,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,E->parallel.world);
 MPI_Allreduce(local+2,out+2,1,MPI_UNSIGNED_LONG_LONG,MPI_BXOR,E->parallel.world);
}
static void write_block(struct All_variables *E,FILE *f,const void *v,size_t n,size_t z) {
 if(fwrite(v,z,n,f)!=n)pices_fail(E,"motion replay write failed");
}
static FILE *open_dump(struct All_variables *E,const char *name) {
 FILE *f=fopen(name,"wb");if(!f)pices_fail(E,"motion replay open failed");return f;
}
static void dump_velocity(struct All_variables *E) {
 char name[1024];int d,n=E->lmesh.nno,step=E->monitor.solution_cycles;
 if(step<81||step>105)return;
 snprintf(name,sizeof(name),"%s.motion_velocity.%d.%d",E->control.data_file,E->parallel.me,step);
 FILE *f=open_dump(E,name);int header[4]={0x504d5631,E->parallel.me,step,n};
 double clock[2]={E->monitor.elapsed_time,E->advection.timestep};
 write_block(E,f,header,4,sizeof(int));write_block(E,f,clock,2,sizeof(double));
 for(d=1;d<=3;d++)write_block(E,f,E->sphere.cap[1].V[d]+1,n,sizeof(float));
 if(fclose(f))pices_fail(E,"motion velocity close failed");
}
void pices_motion_initialize(struct All_variables *E) {
 int p,d,e;char name[1024];if(!pices_motion_enabled(E))return;
 if(E->parallel.nproc>=2097152)pices_fail(E,"motion ID rank range exceeded");
 for(p=1;p<=E->trace.ntracers[1];p++) {
  E->trace.extraq[1][E->pices.slot+1][p]=(double)E->parallel.me*4294967296.0+p;
  E->trace.extraq[1][E->pices.slot+2][p]=0;
 }
 snprintf(name,sizeof(name),"%s.motion_mesh.%d",E->control.data_file,E->parallel.me);
 FILE *f=open_dump(E,name);int head[7]={0x504d4d31,E->parallel.me,E->sphere.capid[1],E->lmesh.nno,E->lmesh.nel,E->lmesh.elx,E->lmesh.elz};
 write_block(E,f,head,7,sizeof(int));
 double geometry[5]={E->sphere.ri,E->sphere.ro,E->trace.box_cushion,E->trace.UV[1][1][0],E->trace.UV[1][2][0]};
 write_block(E,f,geometry,5,sizeof(double));
 for(d=1;d<=3;d++)write_block(E,f,E->sx[1][d]+1,E->lmesh.nno,sizeof(double));
 for(e=1;e<=E->lmesh.nel;e++)write_block(E,f,E->ien[1][e].node+1,8,sizeof(int));
 if(fclose(f))pices_fail(E,"motion mesh close failed");
}
void pices_motion_begin(struct All_variables *E) {
 int p,e,a;char name[1024];if(!pices_motion_enabled(E))return;
 if(!counts) {
  counts=calloc(E->lmesh.nel+1,sizeof(int));selected=calloc(E->lmesh.nel+1,sizeof(int));
  mask=calloc(E->lmesh.nno+1,sizeof(double));
  if(!counts||!selected||!mask)pices_fail(E,"motion mask allocation failed");
  snprintf(name,sizeof(name),"%s.motion.%d",E->control.data_file,E->parallel.me);
  paths=fopen(name,"w");if(!paths)pices_fail(E,"motion log open failed");
  fprintf(paths,"# schema=1 units=length_nd,time_nd,velocity_nd gradient=dv_cartesian_dx_cartesian sparse=start_count_le4_neighbors_only_steps81_105 id=initial_rank_times_2pow32_plus_index\n");
 }
 memset(counts,0,(E->lmesh.nel+1)*sizeof(int));memset(selected,0,(E->lmesh.nel+1)*sizeof(int));
 memset(mask,0,(E->lmesh.nno+1)*sizeof(double));
 for(p=1;p<=E->trace.ntracers[1];p++) {
  e=E->trace.ielement[1][p];if(e<1||e>E->lmesh.nel)pices_fail(E,"motion invalid initial host");counts[e]++;
 }
 for(e=1;e<=E->lmesh.nel;e++)if(counts[e]<=4) {
  selected[e]=1;
  fprintf(paths,"CELL %d %d %d\n",E->monitor.solution_cycles,e,counts[e]);
  for(a=1;a<=8;a++)mask[E->ien[1][e].node[a]]=1;
 }
 if(E->monitor.solution_cycles>=81&&E->monitor.solution_cycles<=105) {
  double *field[NCS];field[1]=mask;E->exchange_node_d(E,field,E->mesh.levmax);
  for(e=1;e<=E->lmesh.nel;e++)for(a=1;a<=8;a++)if(mask[E->ien[1][e].node[a]]>0)selected[e]=1;
 }
 for(p=1;p<=E->trace.ntracers[1];p++)E->trace.extraq[1][E->pices.slot+2][p]=selected[E->trace.ielement[1][p]];
 signatures(E,reference);memset(counters,0,sizeof(counters));memset(maxima,0,sizeof(maxima));
 dump_velocity(E);
}
void pices_motion_particle(struct All_variables *E,int cap,int p,int stage,const double *raw,const double *v) {
 int d;double clip=0,change=0;if(!pices_motion_enabled(E))return;
 for(d=0;d<3;d++) {double diff=raw[d]-E->trace.basicq[cap][d+3][p];clip+=diff*diff;
  diff=v[d]-E->trace.basicq[cap][d+9][p];change+=diff*diff;}
 counters[0]+=(clip>0);if(stage==1)counters[1]+=(change>0);
 maxima[0]=fmax(maxima[0],sqrt(clip));
 if(!E->trace.extraq[cap][E->pices.slot+2][p])return;
 fprintf(paths,"PATH %d %d %.0f %d %.17g",E->monitor.solution_cycles,stage,E->trace.extraq[cap][E->pices.slot+1][p],E->trace.ielement[cap][p],(double)E->advection.timestep);
 for(d=0;d<3;d++)fprintf(paths," %.17g",E->trace.basicq[cap][d+6][p]);
 for(d=0;d<3;d++)fprintf(paths," %.17g",raw[d]);
 for(d=0;d<3;d++)fprintf(paths," %.17g",E->trace.basicq[cap][d+3][p]);
 for(d=0;d<3;d++)fprintf(paths," %.17g",v[d]);
 if(stage==1) {
  double x=E->trace.basicq[cap][6][p],y=E->trace.basicq[cap][7][p],z=E->trace.basicq[cap][8][p],g[10]={0};
  int valid=pices_motion_gradient(E,cap,E->trace.ielement[cap][p],x,y,z,g);
  fprintf(paths," %d",valid);for(d=0;d<10;d++)fprintf(paths," %.17g",g[d]);
  if(valid)maxima[1]=fmax(maxima[1],fabs(g[0]+g[4]+g[8]));else counters[2]++;
 }
 fprintf(paths,"\n");
}
void pices_motion_end(struct All_variables *E,int stage) {
 unsigned long long ids[3],all[3];double maxall[2];int p;
 if(!pices_motion_enabled(E))return;
 signatures(E,ids);
 if(memcmp(ids,reference,sizeof(ids)))pices_fail(E,"motion particle identity signature changed across migration");
 MPI_Allreduce(counters,all,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,E->parallel.world);
 MPI_Allreduce(maxima,maxall,2,MPI_DOUBLE,MPI_MAX,E->parallel.world);
 if(!E->parallel.me)fprintf(E->fp,"PICES_MOTION step=%d stage=%d count=%llu hash_sum=%llu hash_xor=%llu clipped=%llu limited=%llu gradient_unavailable=%llu clip_max_nd=%.17g sampled_div_max_nd=%.17g\n",E->monitor.solution_cycles,stage,ids[0],ids[1],ids[2],all[0],all[1],all[2],maxall[0],maxall[1]);
 for(p=1;p<=E->trace.ntracers[1];p++)if(E->trace.extraq[1][E->pices.slot+2][p])
  fprintf(paths,"HOST %d %d %.0f %d %d\n",E->monitor.solution_cycles,stage,E->trace.extraq[1][E->pices.slot+1][p],E->trace.ielement[1][p],pices_motion_host_valid(E,1,p));
 fflush(paths);fflush(E->fp);memset(counters,0,sizeof(counters));memset(maxima,0,sizeof(maxima));
}
