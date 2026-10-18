/* Test driver only: replay recorded velocities through the unmodified mover. */
#include <stdint.h>
#include <float.h>
#include "pices.h"
static unsigned long long *replay_ids;
static size_t replay_nids;
static int replay_compare_id(const void *a,const void *b) {
 unsigned long long x=*(const unsigned long long *)a,y=*(const unsigned long long *)b;
 return (x>y)-(x<y);
}
static void replay_load_ids(struct All_variables *E) {
 const char *path=getenv("PICES_REPLAY_IDS"); FILE *f; unsigned long long id,previous=0;
 if(!path)return;
 f=fopen(path,"r");if(!f)pices_fail(E,"replay ID list missing");
 while(fscanf(f,"%llu",&id)==1) {
  if(id<=previous || replay_nids>=2000000)pices_fail(E,"replay IDs must be sorted, unique and bounded");
  if(replay_nids%4096==0) {
   void *next=realloc(replay_ids,(replay_nids+4096)*sizeof(*replay_ids));
   if(!next)pices_fail(E,"replay ID allocation failed");replay_ids=next;
  }
  replay_ids[replay_nids++]=id;previous=id;
 }
 if(!feof(f)||!replay_nids)pices_fail(E,"invalid replay ID list");fclose(f);
}
static void replay_output(struct All_variables *E,int step,int factor,double dt,int legacy) {
 char name[2048];FILE *f;int e,p,*count=calloc(E->lmesh.nel+1,sizeof(int));
 unsigned long long hist[10]={0},all[10],n=0,total;
 if(!count)pices_fail(E,"replay census allocation failed");
 if(legacy)snprintf(name,sizeof(name),"%s.replay.%d",E->control.data_file,E->parallel.me);
 else snprintf(name,sizeof(name),"%s.replay.%d.%d",E->control.data_file,E->parallel.me,step);
 f=fopen(name,"w");if(!f)pices_fail(E,"replay output open failed");
 fprintf(f,"# frozen_step=%d substeps=%d total_dt=%.17g\n",step,factor,dt);
 for(p=1;p<=E->trace.ntracers[1];p++) {
  unsigned long long id=(unsigned long long)E->trace.extraq[1][E->pices.slot+1][p];
  e=E->trace.ielement[1][p];if(e<1||e>E->lmesh.nel||!pices_motion_host_valid(E,1,p))pices_fail(E,"replay invalid host");count[e]++;
  if(replay_ids&&!bsearch(&id,replay_ids,replay_nids,sizeof(id),replay_compare_id))continue;
  fprintf(f,"%llu %.17g %.17g %.17g %d\n",id,E->trace.basicq[1][3][p],E->trace.basicq[1][4][p],E->trace.basicq[1][5][p],e);n++;
 }
 if(fclose(f))pices_fail(E,"replay output close failed");
 snprintf(name,sizeof(name),"%s.replay_cells.%d.%d",E->control.data_file,E->parallel.me,step);
 f=fopen(name,"w");if(!f)pices_fail(E,"replay census open failed");
 for(e=1;e<=E->lmesh.nel;e++) {hist[count[e]<9?count[e]:9]++;if(count[e]<=8)fprintf(f,"%d %d\n",e,count[e]);}
 if(fclose(f))pices_fail(E,"replay census close failed");free(count);
 MPI_Allreduce(hist,all,10,MPI_UNSIGNED_LONG_LONG,MPI_SUM,E->parallel.world);
 MPI_Allreduce(&n,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,E->parallel.world);
 if(replay_ids&&total!=replay_nids)pices_fail(E,"replay selected ID count changed");
 if(!E->parallel.me) {
  fprintf(E->fp,"PICES_REPLAY step=%d substeps=%d tracked=%llu histogram=",step,factor,total);
  for(e=0;e<10;e++)fprintf(E->fp,"%s%llu",e?",":"",all[e]);fprintf(E->fp,"\n");fflush(E->fp);
 }
}
static void motion_replay(struct All_variables *E) {
 const char *root=getenv("PICES_REPLAY_ROOT"),*factor=getenv("PICES_REPLAY_SUBSTEPS"),*endenv=getenv("PICES_REPLAY_END_STEP");
 int subdivisions=factor?atoi(factor):0,first=E->monitor.solution_cycles+1,last=endenv?atoi(endenv):first,step,head[4],d;
 double clock[2];char name[2048];FILE *f;
 if(!root||(subdivisions!=1&&subdivisions!=2&&subdivisions!=4)||last<first||last-first>100)pices_fail(E,"invalid replay controls");
 replay_load_ids(E);
 for(step=first;step<=last;step++) {
  snprintf(name,sizeof(name),"%s/DATA/%d/global.motion_velocity.%d.%d",root,E->parallel.me,E->parallel.me,step);
  f=fopen(name,"rb");if(!f)pices_fail(E,"missing frozen velocity");
  if(fread(head,sizeof(int),4,f)!=4||fread(clock,sizeof(double),2,f)!=2||head[0]!=0x504d5631||head[1]!=E->parallel.me||head[2]!=step||head[3]!=E->lmesh.nno)pices_fail(E,"invalid frozen velocity header");
  if(clock[0]!=(double)E->monitor.elapsed_time||!isfinite(clock[1])||clock[1]<=0)pices_fail(E,"replay checkpoint/velocity time mismatch");
  for(d=1;d<=3;d++)if(fread(E->sphere.cap[1].V[d]+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)pices_fail(E,"truncated frozen velocity");
  if(fgetc(f)!=EOF)pices_fail(E,"trailing frozen velocity data");fclose(f);
  E->monitor.solution_cycles=step;E->advection.timestep=clock[1]/subdivisions;
  for(d=0;d<subdivisions;d++)tracer_move_particles(E);
  replay_output(E,step,subdivisions,clock[1],!endenv);
  /* Production clock is float. Match its rounding after each original step. */
  E->monitor.elapsed_time=clock[0]+clock[1];
 }
 free(replay_ids);
}
