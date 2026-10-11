/* Explicit standalone P5 experiments. Disabled defaults leave solver physics
 * untouched. These reconstructed benchmarks are not production initial data. */
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/resource.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "pices.h"
#include "pices_benchmark.h"
static int kind=0,prescribed=0;
static double omega=10000.,length_scale=1.,clocks[2],times[2];
void get_global_shape_fn();
void temperatures_conform_bcs(struct All_variables *);
void p5_parameters(struct All_variables *E) {
 char mode[32];input_string("p5_case",mode,"off",E->parallel.me);
 if(!strcmp(mode,"off"))kind=0;
 else if(!strcmp(mode,"cold"))kind=1;
 else if(!strcmp(mode,"hot"))kind=2;
 else if(!strcmp(mode,"assim"))kind=3;
 else if(!strcmp(mode,"sharp"))kind=4;
 else pices_fail(E,"p5_case must be off/cold/hot/assim/sharp");
 input_boolean("p5_prescribed_velocity",&prescribed,"off",E->parallel.me);
 input_double("p5_omega",&omega,"10000",E->parallel.me);
 input_double("p5_length_scale",&length_scale,"1",E->parallel.me);
 if(!isfinite(omega)||!isfinite(length_scale)||length_scale<=0 || length_scale>4)
  pices_fail(E,"invalid P5 rotation or length scale");
 if(kind) {
  int restart=0,checkpoint=0;
  input_boolean("restart",&restart,"off",E->parallel.me);
  input_boolean("pices_checkpoint",&checkpoint,"off",E->parallel.me);
  if(restart || checkpoint)pices_fail(E,"P5 experiments are fresh runs, not restart/checkpoint modes");
 }
 if(!kind && (prescribed || length_scale!=1))pices_fail(E,"P5 options require an explicit benchmark case");
}
int p5_coupled(void) {return kind && !prescribed;}
double p5_length_scale(void) {return kind?length_scale:1.;}
static double background(struct All_variables *E,double r) {
 return (E->sphere.ro-r)/(E->sphere.ro-E->sphere.ri);
}
static double initial_at(struct All_variables *E,const double x[3],double time) {
 double r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
 if(kind==4)return r<.77?1.:0.;
 double center=kind==1?.82:.68,angle=prescribed?omega*time:0.;
 double dx=x[0]-center*cos(angle),dy=x[1]-center*sin(angle),dz=x[2];
 double taper=sin(M_PI*(r-E->sphere.ri)/(E->sphere.ro-E->sphere.ri));
 return background(E,r)+(kind==1?-.15:.15)*taper*taper*exp(-(dx*dx+dy*dy+dz*dz)/(.12*.12));
}
void p5_initial(struct All_variables *E) {
 int n;double x[3];
 if(!kind)return;
 if(kind==4 && (!E->pices.enabled || !E->pices.consistent_projection || !prescribed || omega!=0))
  pices_fail(E,"sharp stress requires bounded PICES and prescribed zero velocity");
 if(E->control.restart || E->pices.checkpoint || E->sphere.caps_per_proc!=1 || E->sphere.caps!=12)
  pices_fail(E,"P5 benchmarks require a fresh full-sphere run, one cap/rank, checkpoint off");
 if(E->control.ala_pressure_buoyancy)pices_fail(E,"P5 benchmarks require EBA, not strict ALA");
 if(kind==3 && (!E->control.lith_age || !E->control.lith_age_asml))pices_fail(E,"P5 assimilation requires relaxed TA");
 if(prescribed && (E->control.disptn_number!=0 || E->control.Q0!=0 || E->control.lith_age))
  pices_fail(E,"prescribed P5 transport requires zero sources and no TA");
 if(kind!=3)for(n=1;n<=E->lmesh.nno;n++) {
  int j;for(j=0;j<3;j++)x[j]=E->x[1][j+1][n];
  E->T[1][n]=initial_at(E,x,0);
 }
 if(kind==3)for(n=1;n<=E->lmesh.nno;n++) {
  double r=E->sx[1][3][n],f=sin(M_PI*(r-E->sphere.ri)/(E->sphere.ro-E->sphere.ri));
  E->T[1][n]+=.01*f*f*pow(sin(E->sx[1][1][n]),2)*cos(2*E->sx[1][2][n]);
 }
 temperatures_conform_bcs(E);
 fprintf(E->fp,"P5_INIT schema=1 case=%d velocity=%s length_method=min_directional_rms_scaled_v1 length_scale=%.17g omega=%.17g\n",kind,prescribed?"prescribed":"coupled",length_scale,omega);
}
/* Explicit stress input: discontinuous particle temperatures deliberately differ
 * from interpolation of the nodal step. Never used by normal initialization. */
void p5_particle_initial(struct All_variables *E) {
 int p,a,nodes[9];double w[9];
 if(kind!=4)return;
 for(p=1;p<=E->trace.ntracers[1];p++) {
  double radius=0;tracer_temperature_weights(E,1,p,nodes,w);
  for(a=1;a<=8;a++)radius+=w[a]*E->sx[1][3][nodes[a]];
  E->trace.extraq[1][E->pices.slot][p]=radius<.77?1.:0.;
 }
 fprintf(E->fp,"P7_STRESS particle_radial_step=0.77 bounds=0,1\n");
}
int p5_velocity(struct All_variables *E) {
 int n;
 if(!kind || !prescribed)return 0;
 for(n=1;n<=E->lmesh.nno;n++) {
  E->sphere.cap[1].V[1][n]=E->sphere.cap[1].V[3][n]=0;
  E->sphere.cap[1].V[2][n]=omega*E->sx[1][3][n]*sin(E->sx[1][1][n]);
 }
 return 1;
}
void p5_tick(int slot) {if(kind)clocks[slot]=MPI_Wtime();}
void p5_tock(int slot) {if(kind)times[slot]+=MPI_Wtime()-clocks[slot];}
void p5_report(struct All_variables *E) {
 struct Shape_function GN;struct Shape_function_dx dx;struct Shape_function_dA vol;
 struct rusage usage;double rtf[4][9],s[11]={0},sum[11],lo=DBL_MAX,hi=-DBL_MAX,ext[2],peak=0,wrong=0,maxima[5],glob[5],avg[2];
 int e,g,a,n,j;double sign=kind==1?-1.:1.;long long particles=0,allparticles;
 if(!kind)return;
 for(n=1;n<=E->lmesh.nno;n++) {
  double t=E->T[1][n],anom=sign*(t-background(E,E->sx[1][3][n]));
  if(!isfinite(t))pices_fail(E,"P5 nonfinite temperature");
  lo=fmin(lo,t);hi=fmax(hi,t);peak=fmax(peak,anom);wrong=fmax(wrong,-anom);
 }
 for(e=1;e<=E->lmesh.nel;e++) {
  get_global_shape_fn(E,e,&GN,&dx,&vol,0,1,rtf,E->mesh.levmax,1);
  for(g=1;g<=8;g++) {
   double t=0,bg=0,target=0,x[3]={0},v[3]={0},rho=0,cp=0,w=vol.vpt[g]*g_point[g].weight[2],r,anom,pos;
   for(a=1;a<=8;a++) {
    double N=E->N.vpt[GNVINDEX(a,g)];int z;
    n=E->ien[1][e].node[a];z=(n-1)%E->lmesh.noz+1;t+=N*E->T[1][n];
    rho+=N*E->refstate.rho[z];cp+=N*E->refstate.heat_capacity[z];
    {
     double xn[3],rn=0;
     for(j=0;j<3;j++){xn[j]=E->x[1][j+1][n];rn+=xn[j]*xn[j];x[j]+=N*xn[j];v[j]+=N*E->sphere.cap[1].V[j+1][n];}
     bg+=N*background(E,sqrt(rn));
     if(kind!=3)target+=N*initial_at(E,xn,E->monitor.elapsed_time);
    }
   }
   /* Subtract the FE interpolation of the background, so curved-element
    * geometry cannot create a spurious globally distributed anomaly. */
   anom=sign*(t-bg);pos=fmax(anom,0);
   s[0]+=w;s[1]+=w*rho*cp*t;s[2]+=w*fabs(anom);s[3]+=w*anom*anom;
   for(j=0;j<3;j++)s[4+j]+=w*pos*x[j];s[7]+=w*pos;
   for(j=0;j<3;j++){s[8]+=w*v[j]*v[j];s[10]+=w*pos*x[j]*x[j];}
   if(kind!=3){double d=t-target;s[9]+=w*d*d;}
  }
 }
 MPI_Allreduce(s,sum,11,MPI_DOUBLE,MPI_SUM,E->parallel.world);
 MPI_Allreduce(&lo,ext,1,MPI_DOUBLE,MPI_MIN,E->parallel.world);
 MPI_Allreduce(&hi,ext+1,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
 getrusage(RUSAGE_SELF,&usage);
 maxima[0]=peak;maxima[1]=wrong;maxima[2]=times[0];maxima[3]=times[1];maxima[4]=usage.ru_maxrss;
#ifdef __APPLE__
 maxima[4]/=1024.;
#endif
 MPI_Allreduce(maxima,glob,5,MPI_DOUBLE,MPI_MAX,E->parallel.world);
 MPI_Allreduce(times,avg,2,MPI_DOUBLE,MPI_SUM,E->parallel.world);
 if(E->control.tracer)particles=E->trace.ntracers[1];
 MPI_Allreduce(&particles,&allparticles,1,MPI_LONG_LONG_INT,MPI_SUM,E->parallel.world);
 if(E->parallel.me==0) {
  fprintf(E->fp,"P5_METRIC step=%d time=%.17g Tmin=%.17g Tmax=%.17g volume=%.17g energy_proxy=%.17g anomaly_L1=%.17g anomaly_L2=%.17g cx=%.17g cy=%.17g cz=%.17g positive_mass=%.17g anomaly_rms_radius=%.17g peak=%.17g wrong_sign_peak=%.17g vrms=%.17g rotated_initial_L2=%.17g heat_wall_max=%.17g heat_wall_mean=%.17g stokes_wall_max=%.17g stokes_wall_mean=%.17g rss_max_KiB=%.17g particles=%lld\n",
   E->monitor.solution_cycles,(double)E->monitor.elapsed_time,ext[0],ext[1],sum[0],sum[1],sum[2],sqrt(sum[3]/sum[0]),
   sum[7]>0?sum[4]/sum[7]:0,sum[7]>0?sum[5]/sum[7]:0,sum[7]>0?sum[6]/sum[7]:0,sum[7],sum[7]>0?sqrt(fmax(0,sum[10]/sum[7]-(sum[4]*sum[4]+sum[5]*sum[5]+sum[6]*sum[6])/(sum[7]*sum[7]))):0,glob[0],glob[1],sqrt(sum[8]/sum[0]),sqrt(sum[9]/sum[0]),glob[2],avg[0]/E->parallel.nproc,glob[3],avg[1]/E->parallel.nproc,glob[4],allparticles);
  fflush(E->fp);
 }
}
