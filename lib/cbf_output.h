/* Runtime GRD output for Appendix C CBF. Included by Output_gzdir.c so the
 * established Python binding and build target remain unchanged. */
#ifndef CITCOMS_CBF_OUTPUT_H
#define CITCOMS_CBF_OUTPUT_H
#include "cbf_geometry.h"
#include <errno.h>
#include <string.h>
#include <float.h>
#include <sys/stat.h>
#ifdef USE_CBF_NETCDF
#include <netcdf.h>

static int cbf_nc_text(int nc,int var,const char *key,const char *value)
{return nc_put_att_text(nc,var,key,strlen(value),value);}

/* Write a complete temporary classic NetCDF grid; publication is separate. */
static int cbf_write_grd(struct All_variables *E,int step,int top,
                         const double *z,double integral,const char *path)
{
    int nc=-1,xd,yd,xv,yv,zv,dim[2],err=0,i;
    int node_offset=0;
    double x[721],y[361],radius=(top ? E->sphere.ro:E->sphere.ri)*E->data.radius_km*1000;
    double time=E->monitor.elapsed_time;
    double length=E->data.radius_km*1000.0, k0=E->data.k0;
    double deltaT=E->data.ref_temperature;
    double seconds=time*length*length/E->data.kappa0;
    int filtered=E->advection.filter_temperature, assimilated=E->control.lith_age;
#define CBF_NC(call) do {if((err=(call))!=NC_NOERR)goto done;} while(0)
    CBF_NC(nc_create(path,NC_CLOBBER,&nc));
    CBF_NC(nc_def_dim(nc,"x",721,&xd));CBF_NC(nc_def_dim(nc,"y",361,&yd));
    CBF_NC(nc_def_var(nc,"x",NC_DOUBLE,1,&xd,&xv));
    CBF_NC(nc_def_var(nc,"y",NC_DOUBLE,1,&yd,&yv));
    dim[0]=yd;dim[1]=xd;CBF_NC(nc_def_var(nc,"z",NC_FLOAT,2,dim,&zv));
    CBF_NC(cbf_nc_text(nc,xv,"units","degrees_east"));
    CBF_NC(cbf_nc_text(nc,xv,"standard_name","longitude"));
    CBF_NC(cbf_nc_text(nc,yv,"standard_name","latitude"));
    CBF_NC(cbf_nc_text(nc,yv,"units","degrees_north"));
    CBF_NC(cbf_nc_text(nc,zv,"units","W m-2"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"Conventions","CF-1.8"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"method","CBF_GLL_Q1"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"weak_form","physical Galerkin energy; shared thermal kernel"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"mapping_method","ray intersection with Q1 FE face; point sampling, not conservative remapping"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"boundary",top ? "top":"bottom"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"positive_direction",top ? "mantle_to_surface":"core_to_mantle"));
    CBF_NC(cbf_nc_text(nc,NC_GLOBAL,"state",step==0 ?
        "initial diagnostic; initialized Tdot, not an accepted transient step":
        "output T and solver Tdot; assimilation/filter increments are not part of Tdot"));
    CBF_NC(nc_put_att_int(nc,NC_GLOBAL,"node_offset",NC_INT,1,&node_offset));
    CBF_NC(nc_put_att_int(nc,NC_GLOBAL,"step",NC_INT,1,&step));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"time_nd",NC_DOUBLE,1,&time));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"radius_m",NC_DOUBLE,1,&radius));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"time_seconds",NC_DOUBLE,1,&seconds));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"length_scale_m",NC_DOUBLE,1,&length));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"temperature_scale_K",NC_DOUBLE,1,&deltaT));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"conductivity_scale_W_m_K",NC_DOUBLE,1,&k0));
    CBF_NC(nc_put_att_int(nc,NC_GLOBAL,"temperature_filter_enabled",NC_INT,1,&filtered));
    CBF_NC(nc_put_att_int(nc,NC_GLOBAL,"lithosphere_assimilation_enabled",NC_INT,1,&assimilated));
    CBF_NC(nc_put_att_double(nc,NC_GLOBAL,"native_integrated_heat_W",NC_DOUBLE,1,&integral));
    CBF_NC(nc_enddef(nc));
    for(i=0;i<721;++i)x[i]=i*0.5;
    for(i=0;i<361;++i)y[i]=-90+i*0.5;
    CBF_NC(nc_put_var_double(nc,xv,x));CBF_NC(nc_put_var_double(nc,yv,y));
    CBF_NC(nc_put_var_double(nc,zv,z));
done:
    if(nc>=0) {int close_error=nc_close(nc);if(!err)err=close_error;}
    if(err) {fprintf(stderr,"CBF NetCDF: %s: %s\n",path,nc_strerror(err));remove(path);}
#undef CBF_NC
    return err!=NC_NOERR;
}

/* Evaluate the Q1 field on the fixed geographic target grid. A conservative
 * spherical cap bound limits the candidates for each face. Each rank owns
 * elements, not unique boundary nodes; shared-edge hits must agree. */
static int cbf_grid(struct All_variables *E,int top,double **grid,double *heat)
{
    const int nx=721,ny=361,n=nx*ny,lev=E->mesh.levmax,elz=E->lmesh.elz;
    const int active=E->parallel.me_loc[3]==(top ? E->parallel.nprocz-1:0);
    const int side=top ? SIDE_TOP:SIDE_BOTTOM;
    const double pi=3.14159265358979323846,rad=pi/180,spacing=0.5;
    double *low,*high,*glow,*ghigh,*nodal[NCS];
    double xyz[4][3],q[4],dm[4],center[3],norm,dot,angle,maxangle,lon,lat,width;
    double ray[3],value,local_heat=0,global_heat,rr;
    int bad=0,allbad,m,e,a,d,i,j,index,node,lo,hi,left,right,raw_i;
    low=(double *)malloc(n*sizeof(double));high=(double *)malloc(n*sizeof(double));
    glow=(double *)malloc(n*sizeof(double));ghigh=(double *)malloc(n*sizeof(double));
    if(!low||!high||!glow||!ghigh)bad=1;
    for(m=1;m<=E->sphere.caps_per_proc;++m) {
        nodal[m]=(double *)calloc(E->lmesh.nno+2,sizeof(double));
        if(!nodal[m])bad=1;
    }
    MPI_Allreduce(&bad,&allbad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(allbad)goto cleanup;
    for(i=0;i<n;++i) {low[i]=DBL_MAX;high[i]=-DBL_MAX;}
    if(active) for(m=1;m<=E->sphere.caps_per_proc;++m) {
        for(i=1;i<=E->lmesh.nsf;++i) {
            node=top ? E->surf_node[m][i]:E->surf_node[m][i]-E->lmesh.noz+1;
            nodal[m][node]=top ? E->slice.shflux_CBF[m][i]:E->slice.bhflux_CBF[m][i];
        }
        for(e=top ? elz:1;e<=E->lmesh.nel;e+=elz) {
            center[0]=center[1]=center[2]=0;
            for(a=0;a<4;++a) {
                node=E->ien[m][e].node[sidenodes[side][a+1]];q[a]=nodal[m][node];
                for(d=0;d<3;++d) {xyz[a][d]=E->X[lev][m][d+1][node];center[d]+=xyz[a][d];}
            }
            cbf_face_gll_mass(xyz,dm);
            for(a=0;a<4;++a) {local_heat+=dm[a]*q[a];}
            norm=sqrt(center[0]*center[0]+center[1]*center[1]+center[2]*center[2]);
            if(!(norm>0)) {bad=1;continue;}
            for(d=0;d<3;++d)center[d]/=norm;
            maxangle=0;
            for(a=0;a<4;++a) {
                dot=rr=0;for(d=0;d<3;++d) {dot+=xyz[a][d]*center[d];rr+=xyz[a][d]*xyz[a][d];}
                dot/=sqrt(rr);dot=fmax(-1,fmin(1,dot));angle=acos(dot);
                if(angle>maxangle)maxangle=angle;
            }
            maxangle+=1.e-10;
            if(maxangle>=pi/2) {bad=1;continue;}
            lat=asin(fmax(-1,fmin(1,center[2])));lon=atan2(center[1],center[0]);
            if(lon<0)lon+=2*pi;
            lo=(int)floor(((lat-maxangle)/rad+90)/spacing);
            hi=(int)ceil(((lat+maxangle)/rad+90)/spacing);
            lo=lo<0 ? 0:lo;hi=hi>=ny ? ny-1:hi;
            if(fabs(lat)+maxangle>=pi/2) {left=0;right=nx-2;}
            else {
                width=asin(fmin(1,sin(maxangle)/cos(lat)));
                left=(int)floor((lon-width)/rad/spacing);
                right=(int)ceil((lon+width)/rad/spacing);
            }
            for(j=lo;j<=hi;++j)for(raw_i=left;raw_i<=right;++raw_i) {
                i=((raw_i%(nx-1))+(nx-1))%(nx-1);
                lat=(-90+j*spacing)*rad;lon=i*spacing*rad;
                ray[0]=cos(lat)*cos(lon);ray[1]=cos(lat)*sin(lon);ray[2]=sin(lat);
                if(j==0||j==ny-1) {ray[0]=ray[1]=0;ray[2]=j==0 ? -1:1;}
                if(cbf_face_ray_value(xyz,q,ray,&value)) {
                    index=j*nx+i;
                    low[index]=fmin(low[index],value);high[index]=fmax(high[index],value);
                }
            }
        }
    }
    MPI_Allreduce(&bad,&allbad,1,MPI_INT,MPI_MAX,E->parallel.world);
    if(allbad)goto cleanup;
    MPI_Reduce(low,glow,n,MPI_DOUBLE,MPI_MIN,0,E->parallel.world);
    MPI_Reduce(high,ghigh,n,MPI_DOUBLE,MPI_MAX,0,E->parallel.world);
    MPI_Reduce(&local_heat,&global_heat,1,MPI_DOUBLE,MPI_SUM,0,E->parallel.world);
    if(E->parallel.me==0) {
        for(j=0;j<ny;++j) {
            for(i=0;i<nx-1;++i) {
                index=j*nx+i;
                if(glow[index]==DBL_MAX || ghigh[index]==-DBL_MAX ||
                   ghigh[index]-glow[index]>1.e-8*fmax(1,fmax(fabs(glow[index]),fabs(ghigh[index]))))bad=1;
                glow[index]=0.5*(glow[index]+ghigh[index]);
            }
            glow[j*nx+nx-1]=glow[j*nx];
        }
        rr=E->data.radius_km*1000;*heat=global_heat*rr*rr;
    }
    MPI_Bcast(&bad,1,MPI_INT,0,E->parallel.world);allbad=bad;
cleanup:
    free(low);free(high);free(ghigh);
    for(m=1;m<=E->sphere.caps_per_proc;++m)free(nodal[m]);
    if(allbad) {free(glow);*grid=NULL;} else *grid=glow;
    return allbad;
}
#endif

static void cbf_output_grids(struct All_variables *E,int step)
{
    void parallel_process_termination();
#ifdef USE_CBF_NETCDF
    double *grid[2]={NULL,NULL},heat[2]={0,0};
    char tmp[2][256],path[2][256];int top,bad=0;
    int enabled[2]={E->output.cbf_output_bhflux,E->output.cbf_output_shflux};
    for(top=0;top<2;++top)if(enabled[top]) {
        if(cbf_grid(E,top,&grid[top],&heat[top])) {
            if(E->parallel.me==0)fprintf(stderr,"CBF GRD: invalid geometry, uncovered grid or inconsistent shared edge\n");
            parallel_process_termination();
        }
        snprintf(path[top],256,"PostProc/HF_CBF/%s_CBF_%d.grd",top ? "eshf":"cmbhf",step);
        snprintf(tmp[top],256,"%s.tmp",path[top]);
    }
    if(E->parallel.me==0) {
        if(mkdir("PostProc",0775)!=0 && errno!=EEXIST)bad=1;
        if(mkdir("PostProc/HF_CBF",0775)!=0 && errno!=EEXIST)bad=1;
        if(!bad)for(top=0;top<2;++top)if(enabled[top])
            bad|=cbf_write_grd(E,step,top,grid[top],heat[top],tmp[top]);
        if(!bad)for(top=0;top<2;++top)if(enabled[top] && rename(tmp[top],path[top])!=0)bad=1;
        if(bad)for(top=0;top<2;++top)if(enabled[top])remove(tmp[top]);
    }
    MPI_Bcast(&bad,1,MPI_INT,0,E->parallel.world);
    free(grid[0]);free(grid[1]);
    if(bad)parallel_process_termination();
#else
    if(E->parallel.me==0)fprintf(stderr,"CBF GRD requires NetCDF C library; configure with nc-config include/library flags\n");
    parallel_process_termination();
#endif
}
#endif
