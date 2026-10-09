/* P2: legacy payload + strict canonical JSON sidecar + collective manifest.
 * Same executable/partition/physics only. Unknown/missing metadata fails closed.
 * Publication: temporary payload -> temporary metadata -> all-rank rename ->
 * per-rank completion manifest with the same digest of all rank metadata. */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include "global_defs.h"
#include "pices.h"
#include "pices_sha256.h"
#ifndef PICES_SOLVER_COMMIT
#define PICES_SOLVER_COMMIT "unknown"
#endif
#define META_SIZE 4096
static void pathcat(struct All_variables *E,char *out,const char *base,const char *suffix) {
 if(snprintf(out,512,"%s%s",base,suffix)>=512) pices_fail(E,"checkpoint path too long");
}
static void hashfile(struct All_variables *E,const char *path,char hex[65]) {
 FILE *f=fopen(path,"rb");unsigned char b[8192];size_t n;PicesSHA s;
 if(!f)pices_fail(E,"missing checkpoint/input file");pices_sha_init(&s);
 while((n=fread(b,1,sizeof(b),f)))pices_sha_add(&s,b,n);
 if(ferror(f))pices_fail(E,"checkpoint hash read failed");
 if(fclose(f))pices_fail(E,"checkpoint close failed");pices_sha_end(&s,hex);
}
static void readtext(struct All_variables *E,const char *path,char out[META_SIZE]) {
 FILE *f=fopen(path,"rb");size_t n;
 if(!f)pices_fail(E,"missing PICES sidecar/complete manifest");
 n=fread(out,1,META_SIZE-1,f);out[n]=0;
 if(ferror(f) || !feof(f) || strlen(out)!=n)pices_fail(E,"invalid PICES metadata length");
 if(fclose(f))pices_fail(E,"metadata close failed");
}
static void writetext(struct All_variables *E,const char *path,const char *text) {
 FILE *f=fopen(path,"wb");if(!f)pices_fail(E,"cannot write PICES metadata");
 if(fputs(text,f)<0 || fflush(f) || fsync(fileno(f)))pices_fail(E,"metadata write failed");
 if(fclose(f))pices_fail(E,"metadata close failed");
}
/* Older constant-viscosity checkpoints keep their original schema. */
int pices_checkpoint_coupled(struct All_variables *E) {
 return E->viscosity.RHEOL==7 || E->control.vbcs_file ||
        E->control.qvis_mode || E->advection.fixed_timestep<=0;
}
static int schema(struct All_variables *E) {
 return pices_checkpoint_coupled(E)?4:(E->composition.on?3:(E->pices.p4?2:1));
}
/* Identity of the current forcing brackets, independent of relocation. Future
 * forcing coverage belongs to the run's input manifest, not this checkpoint. */
static void forcing_hash(struct All_variables *E,PicesSHA *s) {
 extern float find_age_in_MY(struct All_variables *);
 float age=find_age_in_MY(E);int a=age<0?0:(int)age,b=age<0?0:a+1,k,j;
 const char *prefix[]={E->control.velocity_boundary_file,E->control.lith_age_file,
 E->control.flag_depth_file,E->control.flag_depth_new_file,E->control.tf_file};
 char path[1200],digest[65];
 for(k=0;k<5;k++) {
  if(k==0?!E->control.vbcs_file:!E->control.lith_age)continue;
  for(j=a;j<=((k<3)?b:a);j++) {
   int n=k<2?snprintf(path,sizeof(path),"%s%d.%d",prefix[k],j,E->sphere.capid[1]-1):
               snprintf(path,sizeof(path),"%s%d.xyz",prefix[k],j);
   if(n<0 || n>=sizeof(path))pices_fail(E,"forcing identity path too long");
   hashfile(E,path,digest);pices_sha_add(s,digest,64);
  }
 }
}
static void fingerprint(struct All_variables *E,char out[65]) {
 PicesSHA s;int d;double params[]={E->data.Ttop,E->data.ref_temperature,E->control.Atemp,
 E->control.accuracy,E->control.tole_comp,E->advection.fixed_timestep,
 E->control.reference_conductivity,E->control.kd_upper_prefactor,E->control.kd_lower_prefactor,
 E->pices.no_diffusion,E->mesh.topvbc,E->mesh.botvbc,E->viscosity.RHEOL,
 E->viscosity.TDEPV,E->viscosity.SDEPV,E->pices.max_substeps,E->control.p_iterations,
 E->mesh.toptbc,E->mesh.bottbc,E->control.start_age,E->control.remove_rigid_rotation,
 E->data.scalev,E->data.scalet,E->data.timedir,E->viscosity.MAX,E->viscosity.MIN,
 E->viscosity.max_value,E->viscosity.min_value,E->viscosity.update_allowed,
 E->control.v_steps_low,E->control.augmented_Lagr};
 pices_sha_init(&s);pices_sha_add(&s,params,sizeof(params));
 if(E->pices.consistent_projection){const char method[]="bounded_consistent_v1";pices_sha_add(&s,method,sizeof(method));}
 pices_sha_add(&s,E->viscosity.N0,sizeof(E->viscosity.N0[0])*E->viscosity.num_mat);
 pices_sha_add(&s,E->viscosity.E,sizeof(E->viscosity.E[0])*E->viscosity.num_mat);
 pices_sha_add(&s,E->viscosity.T,sizeof(E->viscosity.T[0])*E->viscosity.num_mat);
 pices_sha_add(&s,E->viscosity.Z,sizeof(E->viscosity.Z[0])*E->viscosity.num_mat);
 for(d=1;d<=3;d++)pices_sha_add(&s,E->x[1][d]+1,E->lmesh.nno*sizeof(E->x[1][d][0]));
 pices_sha_add(&s,E->refstate.rho+1,E->lmesh.noz*sizeof(double));
 pices_sha_add(&s,E->refstate.thermal_expansivity+1,E->lmesh.noz*sizeof(double));
 pices_sha_add(&s,E->refstate.gravity+1,E->lmesh.noz*sizeof(double));
 if(E->pices.p4) {
  int n,i,j;double ta[]={4,E->control.lith_age,E->control.lith_age_asml,E->control.lith_age_time,
   E->control.lith_age_depth,E->control.lith_age_asml_tau_Ma,E->control.lith_age_asml_exp,
   E->control.max_plate_age_Ma,E->refstate.temperature_surface};
  pices_sha_add(&s,ta,sizeof(ta));
  pices_sha_add(&s,E->refstate.Tref+1,E->lmesh.noz*sizeof(double));
  /* Interior TB contains transient targets, not Dirichlet data. */
  for(d=1;d<=3;d++)for(n=1;n<=E->lmesh.nno;n++)
   if(n%E->lmesh.noz==0 || n%E->lmesh.noz==1)
    pices_sha_add(&s,&E->sphere.cap[1].TB[d][n],sizeof(E->sphere.cap[1].TB[d][n]));
  if(E->control.lith_age)for(j=1;j<=E->lmesh.noy;j++)for(i=1;i<=E->lmesh.nox;i++) {
   n=E->lmesh.nxs+i-1+(E->lmesh.nys+j-2)*E->mesh.nox;
   pices_sha_add(&s,&E->age_t[n],sizeof(float));
   pices_sha_add(&s,&E->flag_depth2[n],sizeof(float));
  }
 } else for(d=1;d<=3;d++)pices_sha_add(&s,E->sphere.cap[1].TB[d]+1,E->lmesh.nno*sizeof(E->sphere.cap[1].TB[d][0]));
 if(E->pices.eba) {
  double thermal[]={3,E->data.ks,E->data.radius_km,E->control.Q0,E->control.disptn_number,E->control.surface_temp,
    E->control.eba_formulation,E->control.kT_exponent,E->control.kC_ratio,
    E->control.kd_mantle_thickness_km,E->control.kd_transition_depth_km,
    E->control.kd_upper_linear,E->control.kd_upper_quadratic,
    E->control.kd_lower_linear,E->control.kd_lower_quadratic};
  pices_sha_add(&s,thermal,sizeof(thermal));
  pices_sha_add(&s,E->refstate.heat_capacity+1,E->lmesh.noz*sizeof(double));
  for(d=0;d<PHASE_TRANSITIONS;d++) {
   const struct Phase_transition *p=&E->control.phase[d];
   double phase[]={p->depth,p->density_jump,p->entropy_jump,p->Ra,p->clapeyron,p->transT,p->inv_width};
   pices_sha_add(&s,phase,sizeof(phase));
  }
 } else {
  pices_sha_add(&s,E->pices.K+64,E->lmesh.nel*64*sizeof(double));
  pices_sha_add(&s,E->pices.emass+8,E->lmesh.nel*8*sizeof(double));
 }
 if(E->composition.on) {
  double composition[]={8,E->composition.ibuoy_type,E->composition.ncomp,
   E->composition.ichemical_buoyancy,E->trace.reclassify_flavors,
   E->control.kC_primordial_flavor};
  pices_sha_add(&s,composition,sizeof(composition));
  pices_sha_add(&s,E->composition.buoyancy_ratio,E->composition.ncomp*sizeof(double));
 }
 if(pices_checkpoint_coupled(E)) {
  double coupled[]={4,E->control.qvis_mode,E->control.qvis_cohesion_pa,
   E->control.qvis_friction_angle_rad,E->data.ref_viscosity,E->viscosity.cold_scale,
   E->pices.max_timestep_Ma,E->advection.fine_tune_dt,E->control.NMULTIGRID,
   E->control.NASSEMBLE,E->mesh.levmax,E->control.mg_cycle,
   E->control.v_steps_high,E->control.v_steps_upper,E->control.vbcs_file,
   E->control.precondition,E->control.down_heavy,E->control.up_heavy,
   E->viscosity.SMOOTH,E->viscosity.smooth_cycles,E->viscosity.TDEPV_AVE,
   E->viscosity.EQUIVDD,E->viscosity.equivddopt,E->viscosity.rheol_layers,
   E->viscosity.zlith,E->viscosity.z410,E->viscosity.zlm,E->viscosity.zcmb,
   E->refstate.has_lithostatic_pressure};
  pices_sha_add(&s,coupled,sizeof(coupled));
  if(E->refstate.has_lithostatic_pressure)
   pices_sha_add(&s,E->refstate.lithostatic_pressure_pa+1,E->lmesh.noz*sizeof(double));
  pices_sha_add(&s,E->mat[1]+1,E->lmesh.nel*sizeof(E->mat[1][0]));
  forcing_hash(E,&s);
 }
 pices_sha_end(&s,out);
}
static void metadata(struct All_variables *E,const char *payload,const char *state,char out[META_SIZE]) {
 char sha[65],physics[65],vsha[65];int n;
 hashfile(E,state,vsha);hashfile(E,payload,sha);fingerprint(E,physics);
 n=snprintf(out,META_SIZE,
 "{\n\"magic\":\"CITCOMS_EBA_PICES\",\n\"schema\":%d,\n\"solver_commit\":\"%s\",\n\"phase\":\"accepted\",\n"
 "\"step\":%d,\n\"time\":%.17g,\n\"dt\":%.17g,\n\"total_timesteps\":%d,\n"
 "\"rank\":%d,\n\"mpi_size\":%d,\n\"decomposition\":[%d,%d,%d],\n\"local_mesh\":[%d,%d,%d],\n\"cap_ids\":[%d],\n"
 "\"particles\":%d,\n\"basic_quantities\":%d,\n\"extraq\":[%s{\"name\":\"Tp\",\"slot\":%d}],\n\"flavors\":%d,\n"
 "\"temperature_normalization\":{\"offset_K\":%.17g,\"scale_K\":%.17g},\n"
 "\"interpolation\":\"gnomonic_wedge_radial_v1\",\n\"length\":\"min_directional_rms_cartesian_v1\",\n"
 "\"accepted_velocity_sha256\":\"%s\",\n"
 "\"physics_mesh_sha256\":\"%s\",\n\"checkpoint_sha256\":\"%s\"\n}\n",
 schema(E),PICES_SOLVER_COMMIT,E->monitor.solution_cycles,(double)E->monitor.elapsed_time,(double)E->advection.timestep,E->advection.total_timesteps,
 E->parallel.me,E->parallel.nproc,E->parallel.nprocx,E->parallel.nprocy,E->parallel.nprocz,E->lmesh.nox,E->lmesh.noy,E->lmesh.noz,E->sphere.capid[1],
 E->trace.ntracers[1],E->trace.number_of_basic_quantities,E->trace.nflavors?"{\"name\":\"flavor\",\"slot\":0},":"",E->pices.slot,E->trace.nflavors,
 E->data.Ttop,E->data.ref_temperature,vsha,physics,sha);
 if(n<0 || n>=META_SIZE)pices_fail(E,"metadata overflow");
}
static void manifest(struct All_variables *E,const char *meta,char out[META_SIZE]) {
 PicesSHA s;char digest[65],collective[65],*all;
 hashfile(E,meta,digest);all=malloc(E->parallel.nproc*65);
 if(!all)pices_fail(E,"manifest allocation failed");
 MPI_Allgather(digest,65,MPI_CHAR,all,65,MPI_CHAR,E->parallel.world);
 pices_sha_init(&s);pices_sha_add(&s,all,E->parallel.nproc*65);pices_sha_end(&s,collective);free(all);
 snprintf(out,META_SIZE,"{\"magic\":\"CITCOMS_EBA_PICES_COMPLETE\",\"schema\":1,\"mpi_size\":%d,\"metadata_set_sha256\":\"%s\"}\n",E->parallel.nproc,collective);
}
void pices_checkpoint_publish(struct All_variables *E,const char *temp,const char *final) {
 char mt[512],mf[512],ct[512],cf[512],vt[512],vf[512],text[META_SIZE];
 FILE *f;int d;
 if(strlen(PICES_SOLVER_COMMIT)!=40)pices_fail(E,"checkpoint requires compiled solver commit");
 pathcat(E,mt,final,".pices.json.tmp");pathcat(E,mf,final,".pices.json");
 pathcat(E,ct,final,".pices.manifest.tmp");pathcat(E,cf,final,".pices.manifest");
 /* Remove old completion marker before replacing anything it certified. */
 unlink(cf);
 pathcat(E,vt,final,".pices.state.tmp");pathcat(E,vf,final,".pices.state");
 f=fopen(vt,"wb");if(!f)pices_fail(E,"cannot save accepted velocity");
 for(d=1;d<=3;d++)if(fwrite(E->sphere.cap[1].V[d]+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)pices_fail(E,"velocity write failed");
 if(E->pices.p4) {
  if(fwrite(&E->pices.cbf_valid,sizeof(int),1,f)!=1 ||
     fwrite(E->pices.heat_residual+8,sizeof(double),8*E->lmesh.nel,f)!=(size_t)(8*E->lmesh.nel))
   pices_fail(E,"P4 heat residual write failed");
 }
 if(E->composition.on) {
  int count=E->trace.reclassify_flavors?E->mesh.nox*E->mesh.noy:0;
  if(fwrite(&E->trench_visit_age,sizeof(int),1,f)!=1 || fwrite(&count,sizeof(int),1,f)!=1 ||
     (count && fwrite(E->new_flag_depth+1,sizeof(float),count,f)!=(size_t)count))
   pices_fail(E,"composition history write failed");
 }
 if(pices_checkpoint_coupled(E)) {
  if(fwrite(E->EVI[E->mesh.levmax][1]+1,sizeof(float),8*E->lmesh.nel,f)!=(size_t)(8*E->lmesh.nel) ||
     fwrite(E->VI[E->mesh.levmax][1]+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)
   pices_fail(E,"accepted viscosity write failed");
 }
 if(fflush(f)||fsync(fileno(f))||fclose(f))pices_fail(E,"velocity close failed");
 metadata(E,temp,vt,text);writetext(E,mt,text);
 MPI_Barrier(E->parallel.world);
 if(rename(temp,final) || rename(mt,mf) || rename(vt,vf))pices_fail(E,"checkpoint publish rename failed");
 MPI_Barrier(E->parallel.world);manifest(E,mf,text);writetext(E,ct,text);
 MPI_Barrier(E->parallel.world);
 if(rename(ct,cf))pices_fail(E,"manifest publish failed");
 MPI_Barrier(E->parallel.world);
 fprintf(E->fp,"PICES_CHECKPOINT step=%d phase=accepted schema=%d\n",E->monitor.solution_cycles,schema(E));
}
/* Called before any legacy array allocation/read. All ranks must have a complete
 * set. Payload hash + exact expected size prevents malformed legacy counts. */
void pices_checkpoint_preflight(struct All_variables *E,const char *path) {
 char meta[512],complete[512],state[512],text[META_SIZE],expected[META_SIZE],sha[65],claimed[65],*p;
 FILE *f;struct stat st;int np=-1,h[8],tr[5],sent[4];float times[3];long offset,bytes;
 pathcat(E,meta,path,".pices.json");pathcat(E,complete,path,".pices.manifest");
 manifest(E,meta,expected);readtext(E,complete,text);
 if(strcmp(text,expected))pices_fail(E,"incomplete or mixed checkpoint manifest");
 readtext(E,meta,text);hashfile(E,path,sha);p=strstr(text,"\"checkpoint_sha256\":\"");
 if(!p || sscanf(p,"\"checkpoint_sha256\":\"%64[0-9a-f]\"",claimed)!=1 || strcmp(sha,claimed))pices_fail(E,"checkpoint checksum mismatch");
 pathcat(E,state,path,".pices.state");hashfile(E,state,sha);p=strstr(text,"\"accepted_velocity_sha256\":\"");
 if(!p || sscanf(p,"\"accepted_velocity_sha256\":\"%64[0-9a-f]\"",claimed)!=1 || strcmp(sha,claimed))pices_fail(E,"accepted velocity checksum mismatch");
 if(stat(state,&st) || st.st_size!=(pices_checkpoint_coupled(E)?(8L*E->lmesh.nel+E->lmesh.nno)*sizeof(float):0)+3L*E->lmesh.nno*sizeof(float)+(E->pices.p4?sizeof(int)+8L*E->lmesh.nel*sizeof(double):0)+(E->composition.on?2*sizeof(int)+(E->trace.reclassify_flavors?E->mesh.nox*E->mesh.noy*sizeof(float):0):0))pices_fail(E,"accepted velocity length mismatch");
 p=strstr(text,"\"particles\":");if(!p || sscanf(p,"\"particles\":%d",&np)!=1 || np<0 || np>INT_MAX/128)pices_fail(E,"invalid particle metadata");
 if(sizeof(int)!=4 || sizeof(float)!=4 || sizeof(double)!=8)pices_fail(E,"unsupported checkpoint ABI");
 f=fopen(path,"rb");if(!f)pices_fail(E,"missing checkpoint payload");
 if(fread(h,4,8,f)!=8 || fread(times,4,3,f)!=3)pices_fail(E,"truncated checkpoint header");
 if(h[0]!=E->lmesh.nox || h[1]!=E->lmesh.noy || h[2]!=E->lmesh.noz || h[3]!=E->parallel.nprocx || h[4]!=E->parallel.nprocy || h[5]!=E->parallel.nprocz || h[6]!=1 || h[7]!=E->monitor.solution_cycles_init || !isfinite(times[0]) || !isfinite(times[1]) || !isfinite(times[2]))pices_fail(E,"checkpoint mesh/clock mismatch");
 offset=44L+16+16L*(E->lmesh.nno+1)+16+8+8L*(E->lmesh.npno+1+E->lmesh.neq);
 if(fseek(f,offset,SEEK_SET) || fread(sent,4,4,f)!=4 || fread(tr,4,5,f)!=5)pices_fail(E,"truncated tracer checkpoint");
 if(sent[0]||sent[1]||sent[2]||sent[3]||tr[0]!=12 || tr[1]!=1+(E->trace.nflavors>0) || tr[2]!=E->trace.nflavors || tr[4]!=np)pices_fail(E,"checkpoint tracer schema/count mismatch");
 bytes=offset+16+20+(8L*(6+tr[1])+4)*(np+1L);
 if(E->composition.on) {
  int cs[5];
  if(fseek(f,bytes,SEEK_SET) || fread(cs,4,5,f)!=5 ||
     cs[0] || cs[1] || cs[2] || cs[3] || cs[4]!=E->composition.ncomp)
   pices_fail(E,"checkpoint composition schema mismatch");
  bytes+=20+8L*E->composition.ncomp*(E->lmesh.nel+3L);
 }
 if(fstat(fileno(f),&st) || st.st_size!=bytes)pices_fail(E,"checkpoint length mismatch");
 if(fclose(f))pices_fail(E,"checkpoint close failed");
}
void pices_checkpoint_check_state(struct All_variables *E,const char *path) {
 char meta[512],state[512],actual[META_SIZE],expected[META_SIZE];
 pathcat(E,meta,path,".pices.json");pathcat(E,state,path,".pices.state");readtext(E,meta,actual);metadata(E,path,state,expected);
 if(strcmp(actual,expected)) {
  if(E->parallel.me==0)fprintf(stderr,"PICES_CHECKPOINT expected=%s actual=%s\n",expected,actual);
  pices_fail(E,"checkpoint metadata/physics/normalization/commit mismatch");
 }
}

/* U is the Stokes iterate; nodal V can differ after rigid-rotation removal.
 * Restore the accepted tracer velocity exactly, without a second projection. */
void pices_checkpoint_restore_state(struct All_variables *E,const char *path) {
 char state[512];FILE *f;int d,n;
 pathcat(E,state,path,".pices.state");f=fopen(state,"rb");if(!f)pices_fail(E,"missing accepted velocity");
 for(d=1;d<=3;d++) {
  if(fread(E->sphere.cap[1].V[d]+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)pices_fail(E,"velocity read failed");
  for(n=1;n<=E->lmesh.nno;n++)if(!isfinite(E->sphere.cap[1].V[d][n]))pices_fail(E,"nonfinite restored velocity");
 }
 if(E->pices.p4) {
  E->pices.heat_residual=calloc((E->lmesh.nel+1)*8,sizeof(double));
  if(!E->pices.heat_residual)pices_fail(E,"P4 heat residual allocation failed");
  if(fread(&E->pices.cbf_valid,sizeof(int),1,f)!=1 ||
     fread(E->pices.heat_residual+8,sizeof(double),8*E->lmesh.nel,f)!=(size_t)(8*E->lmesh.nel))
   pices_fail(E,"P4 heat residual read failed");
  if(E->pices.cbf_valid!=(E->monitor.solution_cycles>0))pices_fail(E,"invalid P4 residual lifecycle");
  for(n=8;n<8*(E->lmesh.nel+1);n++)if(!isfinite(E->pices.heat_residual[n]))pices_fail(E,"nonfinite restored heat residual");
 }
 if(E->composition.on) {
  int count,n,expected=E->trace.reclassify_flavors?E->mesh.nox*E->mesh.noy:0;
  if(fread(&E->trench_visit_age,sizeof(int),1,f)!=1 || fread(&count,sizeof(int),1,f)!=1 ||
     count!=expected || (count && fread(E->new_flag_depth+1,sizeof(float),count,f)!=(size_t)count))
   pices_fail(E,"composition history read failed");
  for(n=1;n<=count;n++)if(!isfinite(E->new_flag_depth[n]))pices_fail(E,"invalid composition history");
 }
 if(pices_checkpoint_coupled(E)) {
  float *eta=E->EVI[E->mesh.levmax][1],*vi=E->VI[E->mesh.levmax][1];
  if(fread(eta+1,sizeof(float),8*E->lmesh.nel,f)!=(size_t)(8*E->lmesh.nel) ||
     fread(vi+1,sizeof(float),E->lmesh.nno,f)!=(size_t)E->lmesh.nno)
   pices_fail(E,"accepted viscosity read failed");
  for(n=1;n<=8*E->lmesh.nel;n++)if(!isfinite(eta[n]) || eta[n]<=0)pices_fail(E,"invalid restored integration viscosity");
  for(n=1;n<=E->lmesh.nno;n++)if(!isfinite(vi[n]) || vi[n]<=0)pices_fail(E,"invalid restored nodal viscosity");
 }
 if(fclose(f))pices_fail(E,"velocity close failed");
}
