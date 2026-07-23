/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>
#include <string.h>

void myerror(struct All_variables *,char *);

static float solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max);
static float solve_Ahat_p_fhat_CG(struct All_variables *E,
                                  double **V, double **P, double **F,
                                  double imp, int *steps_max);
static float solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **F,
                                    double imp, int *steps_max);
static float solve_Ahat_p_fhat_ALA_PCG(struct All_variables *E,
                                       double **V, double **P, double **F,
                                       double imp, int *steps_max);
static float solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max);
static double initial_vel_residual(struct All_variables *E,
                                   double **V, double **P, double **F,
                                   double imp);
static double incompressibility_residual(struct All_variables *E,
                                         double **V, double **r);
static void strict_ala_continuity_metrics(struct All_variables *E,
                                          double **V, double **r,
                                          double **div_u, int lev,
                                          double *mass_norm,
                                          double *cancellation_l2);
static void strict_ala_depth_diagnostics(struct All_variables *E,
                                         double **r, double **div_u,
                                         int lev, int iteration);
static void strict_ala_coarse_residual_diagnostics(
    struct All_variables *E, double **r, int lev, int iteration);
static double strict_ala_inner_relative_accuracy(struct All_variables *E,
                                                 double outer_relative_residual);
static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy);
#define ALA_PATCH_MAX_ELEMENTS 8
struct ala_shallow_patch_cache {
    int blocks[NCS];
    unsigned char *size[NCS];
    int *elements[NCS];
    double *chol[NCS];
    double *velocity_node_inverse[NCS];
    int shallow_node_min_z;
    int shallow_node_layers;
};
static void build_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_shallow_patch_cache *cache, int lev);
static void free_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_shallow_patch_cache *cache);
static void apply_ala_pressure_preconditioner(struct All_variables *E,
                                              double **r, double **z,
                                              double **work, int lev,
                                              int iteration,
    struct ala_shallow_patch_cache *patch_cache);
static void apply_ala_coarse_fixed_schur(struct All_variables *E,
                                         double **p, double **Ap,
                                         double **velocity,
                                         double **velocity_rhs,
                                         double **velocity_Ax,
                                         double **velocity_direction,
                                         int lev,
                                         double *residual_reduction);


static int ala_shallow_node_index(struct All_variables *E,
    struct ala_shallow_patch_cache *cache, int node, int lev)
{
    int nz,nx,ny,z,x,y;
    nz=E->lmesh.NOZ[lev];
    nx=E->lmesh.NOX[lev];
    ny=E->lmesh.NOY[lev];
    z=(node-1)%nz+1;
    x=((node-1)/nz)%nx;
    y=(node-1)/(nz*nx);
    if(z<cache->shallow_node_min_z || x<0 || x>=nx || y<0 || y>=ny)
        return -1;
    return (z-cache->shallow_node_min_z)
        +x*cache->shallow_node_layers
        +y*cache->shallow_node_layers*nx;
}


static double ala_shallow_patch_entry(struct All_variables *E,
    struct ala_shallow_patch_cache *cache,
    int e1, int e2, int lev, int m)
{
    int a,b,d1,d2,p1,p2,node1,node2,index;
    double g1,g2,value,*inverse;
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    value=0.0;
    for(a=1;a<=ends;a++) {
        node1=E->IEN[lev][m][e1].node[a];
        for(b=1;b<=ends;b++) {
            node2=E->IEN[lev][m][e2].node[b];
            if(node1!=node2)
                continue;
            index=ala_shallow_node_index(E,cache,node1,lev);
            if(index<0)
                continue;
            inverse=cache->velocity_node_inverse[m]+9*index;
            for(d1=0;d1<dims;d1++) {
                p1=(a-1)*dims+d1;
                g1=E->elt_del[lev][m][e1].g[p1][0]
                    +E->elt_c[lev][m][e1].c[p1][0];
                for(d2=0;d2<dims;d2++) {
                    p2=(b-1)*dims+d2;
                    g2=E->elt_del[lev][m][e2].g[p2][0]
                        +E->elt_c[lev][m][e2].c[p2][0];
                    value += g1*inverse[3*d1+d2]*g2;
                }
            }
        }
    }
    return value;
}


static void build_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_shallow_patch_cache *cache, int lev)
{
    int m,ex,ey,ez,dx,dy,dz,b,i,j,k,n,block_count;
    int elx,ely,elz,shallow_min_ez,radial_groups,e;
    int node,a,row,col,p,node_index,node_count,nz,nx,ny,x,y,z;
    int local_blocks,global_blocks,local_fallback,global_fallback;
    int local_elements,global_elements;
    int local_node_fallback,global_node_fallback;
    double depth_km,sum,pivot,maxdiag,det,shift;
    double local_node_energy[2],global_node_energy[2];
    double *inverse,*element_k;
    double matrix[ALA_PATCH_MAX_ELEMENTS][ALA_PATCH_MAX_ELEMENTS];
    double *L;
    const double pivot_tolerance=1.0e-12;

    memset(cache,0,sizeof(*cache));
    elx=E->lmesh.ELX[lev];
    ely=E->lmesh.ELY[lev];
    elz=E->lmesh.ELZ[lev];
    shallow_min_ez=elz+1;
    for(ez=elz;ez>=1;ez--) {
        depth_km=(1.0-0.5*(E->sx[1][3][ez]+E->sx[1][3][ez+1]))
            *E->data.radius_km;
        if(depth_km>E->control.ala_shallow_patch_depth_km)
            break;
        shallow_min_ez=ez;
    }
    radial_groups=(shallow_min_ez<=elz)
        ? (elz-shallow_min_ez+2)/2 : 0;
    cache->shallow_node_min_z=shallow_min_ez;
    cache->shallow_node_layers=(shallow_min_ez<=elz)
        ? elz-shallow_min_ez+2 : 0;
    block_count=((elx+1)/2)*((ely+1)/2)*radial_groups;
    local_blocks=0;
    local_fallback=0;
    local_elements=0;
    local_node_fallback=0;
    local_node_energy[0]=0.0;
    local_node_energy[1]=0.0;
    nz=E->lmesh.NOZ[lev];
    nx=E->lmesh.NOX[lev];
    ny=E->lmesh.NOY[lev];
    node_count=nx*ny*cache->shallow_node_layers;

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        if(node_count>0) {
            cache->velocity_node_inverse[m]=(double *)calloc(9*node_count,
                                                             sizeof(double));
            if(cache->velocity_node_inverse[m]==NULL)
                myerror(E,"Unable to allocate ALA velocity node blocks");
            for(e=1;e<=E->lmesh.NEL[lev];e++)
                for(a=1;a<=enodes[E->mesh.nsd];a++) {
                    node=E->IEN[lev][m][e].node[a];
                    node_index=ala_shallow_node_index(E,cache,node,lev);
                    if(node_index<0)
                        continue;
                    inverse=cache->velocity_node_inverse[m]+9*node_index;
                    element_k=E->elt_k[lev][m][e].k;
                    for(row=0;row<3;row++)
                        for(col=0;col<3;col++) {
                            p=((a-1)*3+row)*24+(a-1)*3+col;
                            inverse[3*row+col] += element_k[p];
                        }
                }
            for(node_index=0;node_index<node_count;node_index++) {
                inverse=cache->velocity_node_inverse[m]+9*node_index;
                inverse[1]=inverse[3]=0.5*(inverse[1]+inverse[3]);
                inverse[2]=inverse[6]=0.5*(inverse[2]+inverse[6]);
                inverse[5]=inverse[7]=0.5*(inverse[5]+inverse[7]);
                local_node_energy[0] += inverse[0]*inverse[0]
                    +inverse[4]*inverse[4]+inverse[8]*inverse[8];
                local_node_energy[1] += 2.0*(inverse[1]*inverse[1]
                    +inverse[2]*inverse[2]+inverse[5]*inverse[5]);
                maxdiag=max(inverse[0],max(inverse[4],inverse[8]));
                shift=E->control.ala_shallow_patch_regularization
                    *max(maxdiag,1.0e-300);
                inverse[0] += shift;
                inverse[4] += shift;
                inverse[8] += shift;
                det=inverse[0]*(inverse[4]*inverse[8]-inverse[5]*inverse[7])
                    -inverse[1]*(inverse[3]*inverse[8]-inverse[5]*inverse[6])
                    +inverse[2]*(inverse[3]*inverse[7]-inverse[4]*inverse[6]);
                if(!isfinite(det) || det<=1.0e-18*maxdiag*maxdiag*maxdiag) {
                    z=(node_index%cache->shallow_node_layers)
                        +cache->shallow_node_min_z;
                    x=(node_index/cache->shallow_node_layers)%nx;
                    y=node_index/(cache->shallow_node_layers*nx);
                    node=z+x*nz+y*nz*nx;
                    memset(inverse,0,9*sizeof(double));
                    for(row=0;row<3;row++)
                        inverse[3*row+row]=E->BI[lev][m]
                            [E->ID[lev][m][node].doff[row+1]];
                    local_node_fallback++;
                }
                else {
                    matrix[0][0]=(inverse[4]*inverse[8]-inverse[5]*inverse[7])/det;
                    matrix[0][1]=(inverse[2]*inverse[7]-inverse[1]*inverse[8])/det;
                    matrix[0][2]=(inverse[1]*inverse[5]-inverse[2]*inverse[4])/det;
                    matrix[1][0]=(inverse[5]*inverse[6]-inverse[3]*inverse[8])/det;
                    matrix[1][1]=(inverse[0]*inverse[8]-inverse[2]*inverse[6])/det;
                    matrix[1][2]=(inverse[2]*inverse[3]-inverse[0]*inverse[5])/det;
                    matrix[2][0]=(inverse[3]*inverse[7]-inverse[4]*inverse[6])/det;
                    matrix[2][1]=(inverse[1]*inverse[6]-inverse[0]*inverse[7])/det;
                    matrix[2][2]=(inverse[0]*inverse[4]-inverse[1]*inverse[3])/det;
                    for(row=0;row<3;row++)
                        for(col=0;col<3;col++)
                            inverse[3*row+col]=matrix[row][col];
                }
            }
        }
        cache->blocks[m]=block_count;
        if(block_count==0)
            continue;
        cache->size[m]=(unsigned char *)calloc(block_count,
                                                sizeof(unsigned char));
        cache->elements[m]=(int *)calloc(block_count*ALA_PATCH_MAX_ELEMENTS,
                                         sizeof(int));
        cache->chol[m]=(double *)calloc(block_count*ALA_PATCH_MAX_ELEMENTS
                                       *ALA_PATCH_MAX_ELEMENTS,sizeof(double));
        if(cache->size[m]==NULL || cache->elements[m]==NULL ||
           cache->chol[m]==NULL)
            myerror(E,"Unable to allocate ALA shallow-patch cache");
        b=0;
        for(ey=1;ey<=ely;ey+=2)
            for(ex=1;ex<=elx;ex+=2)
                for(ez=elz;ez>=shallow_min_ez;ez-=2) {
                    n=0;
                    for(dy=0;dy<2 && ey+dy<=ely;dy++)
                        for(dx=0;dx<2 && ex+dx<=elx;dx++)
                            for(dz=0;dz<2 && ez-dz>=shallow_min_ez;dz++) {
                                e=(ez-dz)+(ex+dx-1)*elz
                                    +(ey+dy-1)*elz*elx;
                                cache->elements[m][b*ALA_PATCH_MAX_ELEMENTS+n]
                                    =e;
                                n++;
                            }
                    cache->size[m][b]=(unsigned char)n;
                    maxdiag=0.0;
                    for(i=0;i<n;i++)
                        for(j=0;j<n;j++)
                            matrix[i][j]=0.0;
                    for(i=0;i<n;i++) {
                        e=cache->elements[m][b*ALA_PATCH_MAX_ELEMENTS+i];
                        matrix[i][i]=(1.0+E->control.ala_shallow_patch_regularization)
                            *ala_shallow_patch_entry(E,cache,e,e,lev,m);
                        maxdiag=max(maxdiag,matrix[i][i]);
                        for(j=0;j<i;j++) {
                            matrix[i][j]=ala_shallow_patch_entry(
                                E,cache,e,
                                cache->elements[m][b*ALA_PATCH_MAX_ELEMENTS+j],
                                lev,m);
                            matrix[j][i]=matrix[i][j];
                        }
                    }
                    L=cache->chol[m]+b*ALA_PATCH_MAX_ELEMENTS
                        *ALA_PATCH_MAX_ELEMENTS;
                    for(i=0;i<n;i++) {
                        for(j=0;j<=i;j++) {
                            sum=matrix[i][j];
                            for(k=0;k<j;k++)
                                sum -= L[i*ALA_PATCH_MAX_ELEMENTS+k]
                                    *L[j*ALA_PATCH_MAX_ELEMENTS+k];
                            if(i==j) {
                                pivot=sum;
                                if(!isfinite(pivot) ||
                                   pivot<=pivot_tolerance*maxdiag)
                                    break;
                                L[i*ALA_PATCH_MAX_ELEMENTS+j]=sqrt(pivot);
                            }
                            else
                                L[i*ALA_PATCH_MAX_ELEMENTS+j]=sum
                                    /L[j*ALA_PATCH_MAX_ELEMENTS+j];
                        }
                        if(j<=i)
                            break;
                    }
                    if(i<n) {
                        cache->size[m][b]=0;
                        local_fallback++;
                    }
                    else
                        local_blocks++;
                    local_elements += n;
                    b++;
                }
    }
    MPI_Allreduce(&local_blocks,&global_blocks,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_fallback,&global_fallback,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_elements,&global_elements,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_node_fallback,&global_node_fallback,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(local_node_energy,global_node_energy,2,MPI_DOUBLE,MPI_SUM,
                  E->parallel.world);
    if(global_blocks+global_fallback==0)
        myerror(E,"ALA shallow-patch depth selects no global elements");
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA shallow-patch preconditioner depth_km=%e "
                "block=2x2x2 weight=%e regularization=%e "
                "operator=(D+C)Knode^-1(D+C)^T selected_elements=%d "
                "valid_blocks=%d fallback_blocks=%d node_fallbacks=%d "
                "Knode_offdiag_ratio=%e\n",
                E->control.ala_shallow_patch_depth_km,
                E->control.ala_shallow_patch_weight,
                E->control.ala_shallow_patch_regularization,
                global_elements,global_blocks,global_fallback,
                global_node_fallback,sqrt(global_node_energy[1]
                    /max(global_node_energy[0],1.0e-300)));
        fprintf(stderr,"ALA shallow-patch preconditioner depth_km=%e "
                "block=2x2x2 weight=%e regularization=%e "
                "operator=(D+C)Knode^-1(D+C)^T selected_elements=%d "
                "valid_blocks=%d fallback_blocks=%d node_fallbacks=%d "
                "Knode_offdiag_ratio=%e\n",
                E->control.ala_shallow_patch_depth_km,
                E->control.ala_shallow_patch_weight,
                E->control.ala_shallow_patch_regularization,
                global_elements,global_blocks,global_fallback,
                global_node_fallback,sqrt(global_node_energy[1]
                    /max(global_node_energy[0],1.0e-300)));
    }
}


static void free_ala_shallow_patch_cache(struct All_variables *E,
    struct ala_shallow_patch_cache *cache)
{
    int m;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)cache->size[m]);
        free((void *)cache->elements[m]);
        free((void *)cache->chol[m]);
        free((void *)cache->velocity_node_inverse[m]);
    }
    memset(cache,0,sizeof(*cache));
}


static void apply_ala_coarse_fixed_schur(struct All_variables *E,
                                         double **p, double **Ap,
                                         double **velocity,
                                         double **velocity_rhs,
                                         double **velocity_Ax,
                                         double **velocity_direction,
                                         int lev,
                                         double *residual_reduction)
{
    int m,j,k,neq;
    double theta,delta,sigma,rho_cheb,rho_new,residual;
    double local_norm[2],global_norm[2];
    void assemble_grad_rho_p();
    void assemble_div_rho_u();
    void assemble_del2_u();

    neq=E->lmesh.NEQ[lev];
    assemble_grad_rho_p(E,p,velocity_rhs,lev);
    if(strcmp(E->control.ala_two_level_velocity_solver,"diagonal")==0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++)
                velocity[m][j]=E->ALA_velocity_BI[lev][m][j]
                    *velocity_rhs[m][j];
    }
    else {
        theta=0.5*(E->control.ala_two_level_velocity_eigenvalue_max+
                   E->control.ala_two_level_velocity_eigenvalue_min);
        delta=0.5*(E->control.ala_two_level_velocity_eigenvalue_max-
                   E->control.ala_two_level_velocity_eigenvalue_min);
        sigma=theta/delta;
        rho_cheb=1.0/sigma;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++) {
                velocity_direction[m][j]=E->ALA_velocity_BI[lev][m][j]
                    *velocity_rhs[m][j]/theta;
                velocity[m][j]=velocity_direction[m][j];
            }
        for(k=1;k<E->control.ala_two_level_velocity_iterations;k++) {
            assemble_del2_u(E,velocity,velocity_Ax,lev,1);
            rho_new=1.0/(2.0*sigma-rho_cheb);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(j=0;j<neq;j++) {
                    residual=velocity_rhs[m][j]-velocity_Ax[m][j];
                    velocity_direction[m][j]=rho_new*rho_cheb
                        *velocity_direction[m][j]
                        +(2.0*rho_new/delta)*E->ALA_velocity_BI[lev][m][j]
                        *residual;
                    velocity[m][j] += velocity_direction[m][j];
                }
            rho_cheb=rho_new;
        }
    }
    if(residual_reduction!=NULL) {
        assemble_del2_u(E,velocity,velocity_Ax,lev,1);
        local_norm[0]=0.0;
        local_norm[1]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<neq;j++) {
                local_norm[0] += velocity_rhs[m][j]*velocity_rhs[m][j];
                residual=velocity_rhs[m][j]-velocity_Ax[m][j];
                local_norm[1] += residual*residual;
            }
        MPI_Allreduce(local_norm,global_norm,2,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        *residual_reduction=sqrt(global_norm[1]
                                 /max(global_norm[0],1.0e-300));
    }
    assemble_div_rho_u(E,velocity,Ap,lev);
}


static void apply_ala_pressure_preconditioner(struct All_variables *E,
                                              double **r, double **z,
                                              double **work, int lev,
                                              int iteration,
    struct ala_shallow_patch_cache *patch_cache)
{
    int m,j,col,k,e,ex,ey,ez,elz,ncolumns,npno,b,n,i;
    int clev,factor,celx,celz,cnpno,cneq,ce,cx,cy,cz;
    double damping,theta,delta,sigma,rho_cheb,rho_new;
    double velocity_residual_reduction;
    double rhs[ALA_PATCH_MAX_ELEMENTS],solution[ALA_PATCH_MAX_ELEMENTS];
    double sum,weight,*L;
    double local_patch_energy[2],global_patch_energy[2];
    double local_energy[4],global_energy[4];
    double *coarse_rhs[NCS],*coarse_x[NCS],*coarse_residual[NCS];
    double *coarse_Ax[NCS],*coarse_velocity[NCS],*coarse_direction[NCS];
    double *coarse_velocity_rhs[NCS],*coarse_velocity_Ax[NCS];
    double *coarse_velocity_direction[NCS];

    npno=E->lmesh.NPNO[lev];
    if(E->control.ala_radial_line_preconditioner) {
        elz=E->lmesh.ELZ[lev];
        ncolumns=E->lmesh.ELX[lev]*E->lmesh.ELY[lev];
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(col=0;col<ncolumns;col++) {
                if(!E->ALA_BPI_line_valid[lev][m][col+1]) {
                    for(k=0;k<elz;k++) {
                        e=col*elz+k+1;
                        z[m][e]=E->BPI[lev][m][e]*r[m][e];
                    }
                    continue;
                }

                e=col*elz+1;
                work[m][e]=r[m][e];
                for(k=1;k<elz;k++) {
                    e=col*elz+k+1;
                    work[m][e]=r[m][e]-
                        E->ALA_BPI_line_lower[lev][m][e]*work[m][e-1];
                }
                for(k=0;k<elz;k++) {
                    e=col*elz+k+1;
                    z[m][e]=work[m][e]/E->ALA_BPI_line_diag[lev][m][e];
                }
                for(k=elz-2;k>=0;k--) {
                    e=col*elz+k+1;
                    z[m][e]-=E->ALA_BPI_line_lower[lev][m][e+1]*z[m][e+1];
                }
            }
    }
    else {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++)
                z[m][j]=E->BPI[lev][m][j]*r[m][j];
    }

    if(E->control.ala_shallow_patch_preconditioner) {
        weight=E->control.ala_shallow_patch_weight;
        local_patch_energy[0]=0.0;
        local_patch_energy[1]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(b=0;b<patch_cache->blocks[m];b++) {
                n=patch_cache->size[m][b];
                if(n==0)
                    continue;
                L=patch_cache->chol[m]+b*ALA_PATCH_MAX_ELEMENTS
                    *ALA_PATCH_MAX_ELEMENTS;
                for(i=0;i<n;i++) {
                    e=patch_cache->elements[m][b*ALA_PATCH_MAX_ELEMENTS+i];
                    rhs[i]=r[m][e];
                    local_patch_energy[0] += r[m][e]*E->BPI[lev][m][e]
                        *r[m][e];
                    sum=rhs[i];
                    for(j=0;j<i;j++)
                        sum -= L[i*ALA_PATCH_MAX_ELEMENTS+j]*solution[j];
                    solution[i]=sum/L[i*ALA_PATCH_MAX_ELEMENTS+i];
                }
                for(i=n-1;i>=0;i--) {
                    sum=solution[i];
                    for(j=i+1;j<n;j++)
                        sum -= L[j*ALA_PATCH_MAX_ELEMENTS+i]*solution[j];
                    solution[i]=sum/L[i*ALA_PATCH_MAX_ELEMENTS+i];
                }
                for(i=0;i<n;i++) {
                    e=patch_cache->elements[m][b*ALA_PATCH_MAX_ELEMENTS+i];
                    local_patch_energy[1] += r[m][e]*solution[i];
                    z[m][e]=(1.0-weight)*z[m][e]+weight*solution[i];
                }
            }
        if(iteration==0 ||
           iteration%E->control.ala_coarse_residual_interval==0) {
            MPI_Allreduce(local_patch_energy,global_patch_energy,2,MPI_DOUBLE,
                          MPI_SUM,E->parallel.world);
            if(E->parallel.me==0) {
                fprintf(E->fp,"ALA_SHALLOW_PATCH_ENERGY iteration=%d "
                        "diagonal=%e block=%e block_to_diagonal=%e "
                        "blend_weight=%e\n",iteration,
                        global_patch_energy[0],global_patch_energy[1],
                        global_patch_energy[1]
                        /max(global_patch_energy[0],1.0e-300),weight);
                fflush(E->fp);
            }
        }
    }

    if(!E->control.ala_two_level_preconditioner)
        return;

    clev=lev-E->control.ala_two_level_offset;
    factor=1 << E->control.ala_two_level_offset;
    celx=E->lmesh.ELX[clev];
    celz=E->lmesh.ELZ[clev];
    cnpno=E->lmesh.NPNO[clev];
    cneq=E->lmesh.NEQ[clev];
    damping=E->control.ala_two_level_coarse_damping;

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        coarse_rhs[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_x[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_residual[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_Ax[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_velocity[m]=(double *)calloc(cneq+1,sizeof(double));
        coarse_direction[m]=(double *)calloc(cnpno+1,sizeof(double));
        coarse_velocity_rhs[m]=(double *)calloc(cneq+1,sizeof(double));
        coarse_velocity_Ax[m]=(double *)calloc(cneq+1,sizeof(double));
        coarse_velocity_direction[m]=(double *)calloc(cneq+1,sizeof(double));
        if(coarse_rhs[m]==NULL || coarse_x[m]==NULL ||
           coarse_residual[m]==NULL || coarse_Ax[m]==NULL ||
           coarse_velocity[m]==NULL || coarse_direction[m]==NULL ||
           coarse_velocity_rhs[m]==NULL || coarse_velocity_Ax[m]==NULL ||
           coarse_velocity_direction[m]==NULL)
            myerror(E,"Unable to allocate ALA two-level preconditioner");
    }

    /* P^T restriction: sum every factor^3 fine residuals into the
       geometrically coincident coarse pressure cell. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse_rhs[m][ce] += r[m][e];
                }

    /* Both coarse solvers are fixed polynomials of the complete coarse ALA
       operator.  With zero initial state and the BPI similarity transform,
       they define linear symmetric maps for the outer standard PCG. */
    if(strcmp(E->control.ala_two_level_coarse_solver,"jacobi")==0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++)
                coarse_residual[m][ce]=coarse_rhs[m][ce];
        for(k=0;k<E->control.ala_two_level_coarse_iterations;k++) {
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++)
                    coarse_x[m][ce] += damping*E->BPI[clev][m][ce]
                        *coarse_residual[m][ce];
            if(k+1<E->control.ala_two_level_coarse_iterations) {
                apply_ala_coarse_fixed_schur(
                    E,coarse_x,coarse_Ax,coarse_velocity,
                    coarse_velocity_rhs,coarse_velocity_Ax,
                    coarse_velocity_direction,clev,NULL);
                for(m=1;m<=E->sphere.caps_per_proc;m++)
                    for(ce=1;ce<=cnpno;ce++)
                        coarse_residual[m][ce]=coarse_rhs[m][ce]
                            -coarse_Ax[m][ce];
            }
        }
    }
    else {
        theta=0.5*(E->control.ala_two_level_coarse_eigenvalue_max+
                   E->control.ala_two_level_coarse_eigenvalue_min);
        delta=0.5*(E->control.ala_two_level_coarse_eigenvalue_max-
                   E->control.ala_two_level_coarse_eigenvalue_min);
        sigma=theta/delta;
        rho_cheb=1.0/sigma;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(ce=1;ce<=cnpno;ce++) {
                coarse_direction[m][ce]=E->BPI[clev][m][ce]
                    *coarse_rhs[m][ce]/theta;
                coarse_x[m][ce]=coarse_direction[m][ce];
            }
        for(k=1;k<E->control.ala_two_level_coarse_iterations;k++) {
            apply_ala_coarse_fixed_schur(
                E,coarse_x,coarse_Ax,coarse_velocity,
                coarse_velocity_rhs,coarse_velocity_Ax,
                coarse_velocity_direction,clev,NULL);
            rho_new=1.0/(2.0*sigma-rho_cheb);
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(ce=1;ce<=cnpno;ce++) {
                    coarse_residual[m][ce]=coarse_rhs[m][ce]
                        -coarse_Ax[m][ce];
                    coarse_direction[m][ce]=rho_new*rho_cheb
                        *coarse_direction[m][ce]
                        +(2.0*rho_new/delta)*E->BPI[clev][m][ce]
                        *coarse_residual[m][ce];
                    coarse_x[m][ce] += coarse_direction[m][ce];
                }
            rho_cheb=rho_new;
        }
    }

    if(iteration==0 ||
       iteration%E->control.ala_coarse_residual_interval==0) {
        apply_ala_coarse_fixed_schur(
            E,coarse_x,coarse_Ax,coarse_velocity,
            coarse_velocity_rhs,coarse_velocity_Ax,
            coarse_velocity_direction,clev,
            &velocity_residual_reduction);
        local_energy[0]=0.0;
        local_energy[1]=0.0;
        local_energy[2]=0.0;
        local_energy[3]=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(e=1;e<=npno;e++)
                local_energy[0] += r[m][e]*z[m][e];
            for(ce=1;ce<=cnpno;ce++) {
                local_energy[1] += E->control.ala_two_level_coarse_weight
                    *coarse_rhs[m][ce]*coarse_x[m][ce];
                local_energy[2] += coarse_rhs[m][ce]*coarse_rhs[m][ce];
                local_energy[3] += (coarse_rhs[m][ce]-coarse_Ax[m][ce])
                    *(coarse_rhs[m][ce]-coarse_Ax[m][ce]);
            }
        }
        MPI_Allreduce(local_energy,global_energy,4,MPI_DOUBLE,MPI_SUM,
                      E->parallel.world);
        if(E->parallel.me==0) {
            fprintf(E->fp,
                    "ALA_TWO_LEVEL_ENERGY iteration=%d fine=%e coarse=%e "
                    "coarse_to_fine=%e coarse_residual_reduction=%e "
                    "velocity_residual_reduction=%e\n",
                    iteration,global_energy[0],
                    global_energy[1],global_energy[1]
                    /max(global_energy[0],1.0e-300),
                    sqrt(global_energy[3]
                         /max(global_energy[2],1.0e-300)),
                    velocity_residual_reduction);
            fflush(E->fp);
        }
    }

    /* Constant P prolongation, added to the fine diagonal smoother. */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(ey=1;ey<=E->lmesh.ELY[lev];ey++)
            for(ex=1;ex<=E->lmesh.ELX[lev];ex++)
                for(ez=1;ez<=E->lmesh.ELZ[lev];ez++) {
                    e=ez+(ex-1)*E->lmesh.ELZ[lev]
                        +(ey-1)*E->lmesh.ELZ[lev]*E->lmesh.ELX[lev];
                    cx=(ex-1)/factor+1;
                    cy=(ey-1)/factor+1;
                    cz=(ez-1)/factor+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    z[m][e] += E->control.ala_two_level_coarse_weight
                        *coarse_x[m][ce];
                }

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        free((void *)coarse_rhs[m]);
        free((void *)coarse_x[m]);
        free((void *)coarse_residual[m]);
        free((void *)coarse_Ax[m]);
        free((void *)coarse_velocity[m]);
        free((void *)coarse_direction[m]);
        free((void *)coarse_velocity_rhs[m]);
        free((void *)coarse_velocity_Ax[m]);
        free((void *)coarse_velocity_direction[m]);
    }
}


/* Master loop for pressure and (hence) velocity field */

void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{
    void v_from_vector();
    void p_to_nodes();

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    v_from_vector(E);
    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}


static void strict_ala_continuity_metrics(struct All_variables *E,
                                          double **V, double **r,
                                          double **div_u, int lev,
                                          double *mass_norm,
                                          double *cancellation_l2)
{
    int m, e, nel;
    double volume, div_e, c_e, scale_e;
    double local[2], global[2];
    void assemble_div_u();

    assemble_div_u(E, V, div_u, lev);
    local[0] = 0.0;
    local[1] = 0.0;
    nel = E->lmesh.NEL[lev];
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0)
                continue;
            div_e = div_u[m][e];
            c_e = r[m][e] - div_e;
            scale_e = fabs(div_e) + fabs(c_e);
            local[0] += r[m][e] * r[m][e] / volume;
            local[1] += scale_e * scale_e / volume;
        }

    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    *mass_norm = sqrt(max(global[0], 0.0));
    *cancellation_l2 = sqrt(global[0] / max(global[1], 1.0e-300));
}


static void strict_ala_depth_diagnostics(struct All_variables *E,
                                         double **r, double **div_u,
                                         int lev, int iteration)
{
    int m, e, b, bins, nel, local_ez, global_ez, global_elz;
    int fields;
    double volume, div_e, c_e, scale_e, depth_km, r2;
    double cancellation, correlation, d_over_c, residual_fraction;
    double *local, *global;
    double local_summary[5], global_summary[5];

    if(!E->control.ala_depth_diagnostics)
        return;
    if(iteration != 0 &&
       (iteration % E->control.ala_depth_diagnostic_interval) != 0)
        return;

    global_elz = E->mesh.ELZ[lev];
    bins = min(E->control.ala_depth_diagnostic_bins, global_elz);
    fields = 7;
    local = (double *)calloc(bins * fields, sizeof(double));
    global = (double *)calloc(bins * fields, sizeof(double));
    if(local == NULL || global == NULL)
        myerror(E, "Unable to allocate ALA depth diagnostics");

    for(b=0; b<5; b++)
        local_summary[b] = 0.0;

    nel = E->lmesh.NEL[lev];
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=nel; e++) {
            volume = E->eco[m][e].area;
            if(volume <= 0.0)
                continue;
            local_ez = (e - 1) % E->lmesh.ELZ[lev] + 1;
            global_ez = E->lmesh.EZS[lev] + local_ez;
            b = ((global_ez - 1) * bins) / global_elz;
            b = min(max(b, 0), bins - 1);

            depth_km = (1.0 - 0.5 *
                        (E->sx[m][3][local_ez] +
                         E->sx[m][3][local_ez+1]))
                * E->data.radius_km;
            div_e = div_u[m][e];
            c_e = r[m][e] - div_e;
            scale_e = fabs(div_e) + fabs(c_e);
            r2 = r[m][e] * r[m][e] / volume;

            local[b*fields] += div_e * div_e / volume;
            local[b*fields+1] += c_e * c_e / volume;
            local[b*fields+2] += r2;
            local[b*fields+3] += div_e * c_e / volume;
            local[b*fields+4] += scale_e * scale_e / volume;
            local[b*fields+5] += volume;
            local[b*fields+6] += depth_km * volume;

            local_summary[0] += r2;
            local_summary[1] += depth_km * r2;
            if(depth_km <= 200.0)
                local_summary[2] += r2;
            if(depth_km <= 410.0)
                local_summary[3] += r2;
            if(depth_km <= 660.0)
                local_summary[4] += r2;
        }

    MPI_Allreduce(local, global, bins * fields, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(local_summary, global_summary, 5, MPI_DOUBLE, MPI_SUM,
                  E->parallel.world);

    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA_DEPTH_CONTINUITY iteration=%d solver=%s bins=%d "
                "residual_mean_depth_km=%e top200_fraction=%e "
                "top410_fraction=%e top660_fraction=%e\n",
                iteration, E->control.uzawa, bins,
                global_summary[1] / max(global_summary[0], 1.0e-300),
                global_summary[2] / max(global_summary[0], 1.0e-300),
                global_summary[3] / max(global_summary[0], 1.0e-300),
                global_summary[4] / max(global_summary[0], 1.0e-300));
        fprintf(E->fp,
                "ALA_DEPTH_BIN bin depth_km residual_fraction cancellation "
                "D_over_C D_C_correlation\n");
        for(b=bins-1; b>=0; b--) {
            residual_fraction = global[b*fields+2]
                / max(global_summary[0], 1.0e-300);
            cancellation = sqrt(global[b*fields+2]
                                / max(global[b*fields+4], 1.0e-300));
            d_over_c = sqrt(global[b*fields]
                            / max(global[b*fields+1], 1.0e-300));
            correlation = global[b*fields+3]
                / sqrt(max(global[b*fields] * global[b*fields+1],
                           1.0e-300));
            correlation = min(max(correlation, -1.0), 1.0);
            depth_km = global[b*fields+6]
                / max(global[b*fields+5], 1.0e-300);
            fprintf(E->fp,
                    "ALA_DEPTH_BIN %d %e %e %e %e %e\n",
                    b, depth_km, residual_fraction, cancellation,
                    d_over_c, correlation);
        }
        fflush(E->fp);
    }

    free((void *)local);
    free((void *)global);
}


static void strict_ala_coarse_residual_diagnostics(
    struct All_variables *E, double **r, int lev, int iteration)
{
    int m, e, ex, ey, ez, offset, factor;
    int elx, ely, elz, parent_elx, parent_ely, parent_elz, parent_nel;
    int parent_ex, parent_ey, parent_e;
    double volume, fine_energy, coarse_energy, coarse_fraction;
    double local_fine_energy, global_fine_energy;
    double local_coarse_energy, global_coarse_energy;
    double *parent_r, *parent_volume;

    if(!E->control.ala_coarse_residual_diagnostics)
        return;
    if(iteration != 0 &&
       (iteration % E->control.ala_coarse_residual_interval) != 0)
        return;

    elx = E->lmesh.ELX[lev];
    ely = E->lmesh.ELY[lev];
    elz = E->lmesh.ELZ[lev];
    local_fine_energy = 0.0;
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.NEL[lev]; e++) {
            volume = E->eco[m][e].area;
            if(volume > 0.0)
                local_fine_energy += r[m][e] * r[m][e] / volume;
        }
    MPI_Allreduce(&local_fine_energy, &global_fine_energy, 1, MPI_DOUBLE,
                  MPI_SUM, E->parallel.world);

    if(E->parallel.me == 0)
        fprintf(E->fp,
                "ALA_COARSE_RESIDUAL iteration=%d solver=%s fine_energy=%e\n",
                iteration, E->control.uzawa, global_fine_energy);

    factor = 1;
    for(offset=1; offset<=E->control.ala_coarse_residual_levels; offset++) {
        factor *= 2;
        if(elx % factor || ely % factor || elz % factor)
            myerror(E, "ALA coarse residual levels require locally divisible "
                    "element dimensions");

        parent_elx = elx / factor;
        parent_ely = ely / factor;
        parent_elz = elz / factor;
        parent_nel = parent_elx * parent_ely * parent_elz;
        parent_r = (double *)calloc(parent_nel, sizeof(double));
        parent_volume = (double *)calloc(parent_nel, sizeof(double));
        if(parent_r == NULL || parent_volume == NULL)
            myerror(E, "Unable to allocate ALA coarse residual diagnostics");

        local_coarse_energy = 0.0;
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            memset(parent_r, 0, parent_nel * sizeof(double));
            memset(parent_volume, 0, parent_nel * sizeof(double));
            for(ey=1; ey<=ely; ey++)
                for(ex=1; ex<=elx; ex++)
                    for(ez=1; ez<=elz; ez++) {
                        e = ez + (ex-1)*elz + (ey-1)*elz*elx;
                        volume = E->eco[m][e].area;
                        if(volume <= 0.0)
                            continue;
                        parent_ex = (ex-1) / factor;
                        parent_ey = (ey-1) / factor;
                        parent_e = (ez-1) / factor
                            + parent_ex*parent_elz
                            + parent_ey*parent_elz*parent_elx;
                        parent_r[parent_e] += r[m][e];
                        parent_volume[parent_e] += volume;
                    }
            for(parent_e=0; parent_e<parent_nel; parent_e++)
                if(parent_volume[parent_e] > 0.0)
                    local_coarse_energy += parent_r[parent_e]
                        * parent_r[parent_e] / parent_volume[parent_e];
        }
        MPI_Allreduce(&local_coarse_energy, &global_coarse_energy, 1,
                      MPI_DOUBLE, MPI_SUM, E->parallel.world);

        fine_energy = max(global_fine_energy, 0.0);
        coarse_energy = max(global_coarse_energy, 0.0);
        coarse_fraction = (fine_energy > 0.0)
            ? coarse_energy / fine_energy : 0.0;
        if(coarse_fraction < -1.0e-12 || coarse_fraction > 1.0+1.0e-12)
            myerror(E, "ALA coarse residual fraction is outside [0,1]");
        coarse_fraction = min(max(coarse_fraction, 0.0), 1.0);
        if(E->parallel.me == 0)
            fprintf(E->fp,
                    "ALA_COARSE_LEVEL offset=%d coarse_fraction=%e "
                    "coarse_energy=%e\n",
                    offset, coarse_fraction, coarse_energy);

        free((void *)parent_r);
        free((void *)parent_volume);
    }
    if(E->parallel.me == 0)
        fflush(E->fp);
}

void solve_constrained_flow_iterative_pseudo_surf(E)
     struct All_variables *E;

{
    void v_from_vector_pseudo_surf();
    void p_to_nodes();

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    v_from_vector_pseudo_surf(E);
    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}


/* ========================================================================= */

static void print_convergence_progress(struct All_variables *E,
                                       int count, double time0,
                                       double sq_vdotv,
                                       double dvelocity, double dpressure)
{
    double CPU_time0();

    fprintf(E->fp, "AhatP (%03d) after %6.2f s v=%.3e  div/v=%.3e "
            "dv/v=%.3e and dp/p=%.3e for step %d\n",
            count, CPU_time0()-time0, sq_vdotv,E->monitor.incompressibility,
            dvelocity, dpressure, E->monitor.solution_cycles);
    fprintf(stderr, "AhatP (%03d) after %6.2f s v=%.3e div/v=%.3e "
            "dv/v=%.3e and dp/p=%.3e for step %d\n",
            count, CPU_time0()-time0, sq_vdotv, E->monitor.incompressibility,
            dvelocity, dpressure, E->monitor.solution_cycles);

    return;
}


static double strict_ala_inner_relative_accuracy(
    struct All_variables *E, double outer_relative_residual)
{
    double candidate, final_accuracy, relative_accuracy;

    /* The user cap may deliberately request a tighter K solve than the
     * conventional 0.1*outer-tolerance forcing target. */
    final_accuracy = E->control.ala_inner_accuracy_max;
    if(E->control.tole_comp > 0.0)
        final_accuracy = min(final_accuracy,
                             0.1 * E->control.tole_comp);
    candidate = E->control.ala_inner_accuracy_factor
        * max(outer_relative_residual, 0.0);

    if(candidate <= final_accuracy)
        relative_accuracy = final_accuracy;
    else
        /* Quantize by decades so the inexact Schur operator only changes at
         * explicit Krylov restart points. */
        relative_accuracy = pow(10.0, ceil(log10(candidate) - 1.0e-12));

    relative_accuracy = min(relative_accuracy,
                            E->control.ala_inner_accuracy_max);
    return max(relative_accuracy, final_accuracy);
}


static double strict_ala_inner_accuracy(struct All_variables *E,
                                        double **F, int lev,
                                        double relative_accuracy)
{
    double global_vdot();
    double rhs_norm;

    rhs_norm = sqrt(global_vdot(E, F, F, lev) / E->mesh.neq);

    return max(relative_accuracy * rhs_norm, 1.0e-14);
}



static float solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max)
{
    float residual;

    if(E->control.inv_gruneisen == 0)
        residual = solve_Ahat_p_fhat_CG(E, V, P, F, imp, steps_max);
    else {
        if(strcmp(E->control.uzawa, "cg") == 0)
            residual = solve_Ahat_p_fhat_iterCG(E, V, P, F, imp, steps_max);
        else if(strcmp(E->control.uzawa, "ala_cg") == 0)
            residual = solve_Ahat_p_fhat_ALA_PCG(E, V, P, F, imp,
                                                steps_max);
        else if(strcmp(E->control.uzawa, "bicg") == 0)
            residual = solve_Ahat_p_fhat_BiCG(E, V, P, F, imp, steps_max);
        else
            myerror(E, "Error: unknown Uzawa iteration\n");
    }

    return(residual);
}


/* Solve incompressible Stokes flow using
 * conjugate gradient (CG) iterations
 */

static float solve_Ahat_p_fhat_CG(struct All_variables *E,
                                  double **V, double **P, double **FF,
                                  double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;
    int gnpno, gneq;

    double *r1[NCS], *r2[NCS], *z1[NCS], *s1[NCS], *s2[NCS], *F[NCS];
    double *shuffle[NCS];
    double alpha, delta, r0dotz0, r1dotz1,sq_vdotv;
    double residual, v_res;

    double global_vdot(), global_pdot();

    double time0, CPU_time0();
    float dpressure, dvelocity;

    void assemble_c_u();
    void assemble_div_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();

    gnpno = E->mesh.npno;
    gneq = E->mesh.neq;
    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m] = (double *)malloc(neq*sizeof(double));
        r1[m] = (double *)malloc((npno+1)*sizeof(double));
        r2[m] = (double *)malloc((npno+1)*sizeof(double));
        z1[m] = (double *)malloc((npno+1)*sizeof(double));
        s1[m] = (double *)malloc((npno+1)*sizeof(double));
        s2[m] = (double *)malloc((npno+1)*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;


    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* calculate the contribution of compressibility in the continuity eqn */
    if(E->control.inv_gruneisen != 0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++)
                r2[m][j] = 0.0;

        assemble_c_u(E, V, r2, lev);
    }


    /* calculate the initial velocity residual */
    v_res = initial_vel_residual(E, V, P, F, imp);


    /* initial residual r1 = div(V) */
    assemble_div_u(E, V, r1, lev);


    /* add the contribution of compressibility to the initial residual */
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=1;j<=npno;j++) {
                r1[m][j] += r2[m][j];
            }

    residual = incompressibility_residual(E, V, r1);


    sq_vdotv = sqrt(E->monitor.vdotv);

    /* pressure and velocity corrections */
    dpressure = 1.0;
    dvelocity = 1.0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
    }


    r0dotz0 = 0;

    while( (count < *steps_max) &&
           (E->monitor.incompressibility >= E->control.tole_comp) &&
           (dpressure >= imp) && (dvelocity >= imp) )  {


        /* preconditioner BPI ~= inv(K), z1 = BPI*r1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                z1[m][j] = E->BPI[lev][m][j] * r1[m][j];


        /* r1dotz1 = <r1, z1> */
        r1dotz1 = global_pdot(E, r1, z1, lev);
        assert(r1dotz1 != 0.0  /* Division by zero in head of incompressibility iteration */);


        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j];
        else {
            /* s2 = z1 + s1 * <r1,z1>/<r0,z0> */
            delta = r1dotz1 / r0dotz0;
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    s2[m][j] = z1[m][j] + delta * s1[m][j];
        }


        /* solve K*u1 = grad(s2) for u1 */
        assemble_grad_p(E, s2, F, lev);
        valid = solve_del2_u(E, E->u1, F, imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
	    fprintf(stderr,"Solver not converging! \n");
	    fflush(stderr);
	    exit(18);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* The outer CG path follows CitcomS' split compressible strategy:
         * C*u is held on the outer-loop right hand side. */
        assemble_div_u(E, E->u1, F, lev);


        /* alpha = <r1, z1> / <s2, F> */
        if(valid)
            /* alpha defined this way is the same as R&W */
            alpha = r1dotz1 / global_pdot(E, s2, F, lev);
        else
            alpha = 0.0;


        /* r2 = r1 - alpha * div(u1) */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = r1[m][j] - alpha * F[m][j];


        /* P = P + alpha * s2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += alpha * s2[m][j];


        /* V = V - alpha * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= alpha * E->u1[m][j];


        /* compute velocity and incompressibility residual */
        assemble_div_u(E, V, F, lev);
        incompressibility_residual(E, V, F);

        /* compute velocity and pressure corrections */
        dpressure = alpha * sqrt(global_pdot(E, s2, s2, lev)
                                 / (1.0e-32 + global_pdot(E, P, P, lev)));
        dvelocity = alpha * sqrt(global_vdot(E, E->u1, E->u1, lev)
                                 / (1.0e-32 + E->monitor.vdotv));

        count++;

	sq_vdotv = sqrt(E->monitor.vdotv);

        if (E->control.print_convergence && E->parallel.me==0)  {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
        }


        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = s1[m];
            s1[m] = s2[m];
            s2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, z0> = <r1, z1> */
        r0dotz0 = r1dotz1;

    } /* end loop for conjugate gradient */

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free((void *) F[m]);
        free((void *) r1[m]);
        free((void *) r2[m]);
        free((void *) z1[m]);
        free((void *) s1[m]);
        free((void *) s2[m]);
    }

    *steps_max=count;

    return(residual);
}

/* Solve strict ALA pressure correction using preconditioned conjugate
 * gradients on S=(D+C) K^-1 (D+C)^T.  Unlike the legacy compressible CG
 * path, this applies the complete continuity operator and its transpose.
 */

static float solve_Ahat_p_fhat_ALA_PCG(struct All_variables *E,
                                       double **V, double **P, double **FF,
                                       double imp, int *steps_max)
{
    void assemble_div_rho_u();
    void assemble_grad_rho_p();
    void strip_bcs_from_residual();
    void parallel_process_termination();
    int solve_del2_u();

    double global_vdot(), global_pdot();
    double CPU_time0();

    int npno, neq, lev, coarse_lev, coarse_factor;
    int m, j, count, valid;
    int restart_search;
    int hybrid_consecutive_count, hybrid_converged;
    int nonpositive_curvature_count;
    int local_invalid_bpi, global_invalid_bpi;
    int local_invalid_coarse, global_invalid_coarse;

    double alpha, beta, rho, rho_old, curvature;
    double min_curvature, sq_vdotv;
    double residual, dpressure, dvelocity;
    double initial_rnorm, rnorm, relative_residual;
    double recursive_rnorm, recursive_relative_residual, drift_ratio;
    double initial_mass_norm, mass_norm, mass_relative_residual;
    double cancellation_l2;
    double inner_accuracy, inner_relative_accuracy;
    double local_bpi_min, local_bpi_max;
    double global_bpi_min, global_bpi_max;
    double time0;

    double *F[NCS];
    double *r[NCS], *z[NCS], *p[NCS], *q[NCS];
    double *explicit_r[NCS], *div_u[NCS];
    double *preconditioner_work[NCS];
    struct ala_shallow_patch_cache patch_cache;

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;
    memset(&patch_cache,0,sizeof(patch_cache));

    if(E->control.ala_two_level_preconditioner) {
        coarse_lev=lev-E->control.ala_two_level_offset;
        if(coarse_lev<E->mesh.levmin)
            myerror(E,"ALA two-level offset is below the available mesh hierarchy");
        coarse_factor=1 << E->control.ala_two_level_offset;
        if(E->lmesh.ELX[lev]!=coarse_factor*E->lmesh.ELX[coarse_lev] ||
           E->lmesh.ELY[lev]!=coarse_factor*E->lmesh.ELY[coarse_lev] ||
           E->lmesh.ELZ[lev]!=coarse_factor*E->lmesh.ELZ[coarse_lev])
            myerror(E,"ALA two-level pressure transfer does not match mesh hierarchy");
    }

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        F[m] = (double *)malloc(neq*sizeof(double));
        r[m] = (double *)malloc((npno+1)*sizeof(double));
        z[m] = (double *)malloc((npno+1)*sizeof(double));
        p[m] = (double *)malloc((npno+1)*sizeof(double));
        q[m] = (double *)malloc((npno+1)*sizeof(double));
        explicit_r[m] = (double *)malloc((npno+1)*sizeof(double));
        div_u[m] = (double *)malloc((npno+1)*sizeof(double));
        preconditioner_work[m] =
            (double *)malloc((npno+1)*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;

    local_bpi_min = 1.0e300;
    local_bpi_max = 0.0;
    local_invalid_bpi = 0;
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=1; j<=npno; j++) {
            if(!isfinite(E->BPI[lev][m][j]) || E->BPI[lev][m][j] <= 0.0)
                local_invalid_bpi++;
            else {
                local_bpi_min = min(local_bpi_min, E->BPI[lev][m][j]);
                local_bpi_max = max(local_bpi_max, E->BPI[lev][m][j]);
            }
        }
    MPI_Allreduce(&local_bpi_min, &global_bpi_min, 1, MPI_DOUBLE, MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_bpi_max, &global_bpi_max, 1, MPI_DOUBLE, MPI_MAX,
                  E->parallel.world);
    MPI_Allreduce(&local_invalid_bpi, &global_invalid_bpi, 1, MPI_INT, MPI_SUM,
                  E->parallel.world);
    local_invalid_coarse=0;
    if(E->control.ala_two_level_preconditioner)
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(j=1;j<=E->lmesh.NPNO[coarse_lev];j++)
                if(!isfinite(E->BPI[coarse_lev][m][j]) ||
                   E->BPI[coarse_lev][m][j]<=0.0)
                    local_invalid_coarse++;
            for(j=0;j<E->lmesh.NEQ[coarse_lev];j++)
                if(!isfinite(E->ALA_velocity_BI[coarse_lev][m][j]) ||
                   E->ALA_velocity_BI[coarse_lev][m][j]<=0.0)
                    local_invalid_coarse++;
        }
    MPI_Allreduce(&local_invalid_coarse,&global_invalid_coarse,1,MPI_INT,
                  MPI_SUM,E->parallel.world);
    if(global_invalid_bpi || global_invalid_coarse)
        parallel_process_termination();
    if(E->control.ala_shallow_patch_preconditioner)
        build_ala_shallow_patch_cache(E,&patch_cache,lev);
    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA PCG pressure preconditioner = %s mode=%s "
                "BPI_range=(%e,%e) invalid=%d restart_interval=%d\n",
                E->control.precondition ? "on" : "off",
                E->control.ala_two_level_preconditioner
                    ? "diagonal_plus_coarse"
                    : (E->control.ala_shallow_patch_preconditioner
                       ? "shallow_patch"
                       : (E->control.ala_radial_line_preconditioner
                          ? "radial_line" : "diagonal")),
                global_bpi_min, global_bpi_max, global_invalid_bpi,
                E->control.ala_pcg_restart_interval);
        fprintf(stderr,
                "ALA PCG pressure preconditioner = %s mode=%s "
                "BPI_range=(%e,%e) invalid=%d restart_interval=%d\n",
                E->control.precondition ? "on" : "off",
                E->control.ala_two_level_preconditioner
                    ? "diagonal_plus_coarse"
                    : (E->control.ala_shallow_patch_preconditioner
                       ? "shallow_patch"
                       : (E->control.ala_radial_line_preconditioner
                          ? "radial_line" : "diagonal")),
                global_bpi_min, global_bpi_max, global_invalid_bpi,
                E->control.ala_pcg_restart_interval);
        if(E->control.ala_two_level_preconditioner) {
            fprintf(E->fp,
                    "ALA two-level pressure correction offset=%d level=%d "
                    "coarse_solver=%s coarse_iterations=%d damping=%e "
                    "eigen_interval=(%e,%e) coarse_weight=%e "
                    "velocity_solver=%s velocity_iterations=%d "
                    "velocity_eigen_interval=(%e,%e) "
                    "invalid_diagonals=%d "
                    "operator="
                    "(D+C)Kfixed^-1(D+C)^T transfer=Pt/P\n",
                    E->control.ala_two_level_offset,coarse_lev,
                    E->control.ala_two_level_coarse_solver,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_coarse_damping,
                    E->control.ala_two_level_coarse_eigenvalue_min,
                    E->control.ala_two_level_coarse_eigenvalue_max,
                    E->control.ala_two_level_coarse_weight,
                    E->control.ala_two_level_velocity_solver,
                    E->control.ala_two_level_velocity_iterations,
                    E->control.ala_two_level_velocity_eigenvalue_min,
                    E->control.ala_two_level_velocity_eigenvalue_max,
                    global_invalid_coarse);
            fprintf(stderr,
                    "ALA two-level pressure correction offset=%d level=%d "
                    "coarse_solver=%s coarse_iterations=%d damping=%e "
                    "eigen_interval=(%e,%e) coarse_weight=%e "
                    "velocity_solver=%s velocity_iterations=%d "
                    "velocity_eigen_interval=(%e,%e) "
                    "invalid_diagonals=%d "
                    "operator="
                    "(D+C)Kfixed^-1(D+C)^T transfer=Pt/P\n",
                    E->control.ala_two_level_offset,coarse_lev,
                    E->control.ala_two_level_coarse_solver,
                    E->control.ala_two_level_coarse_iterations,
                    E->control.ala_two_level_coarse_damping,
                    E->control.ala_two_level_coarse_eigenvalue_min,
                    E->control.ala_two_level_coarse_eigenvalue_max,
                    E->control.ala_two_level_coarse_weight,
                    E->control.ala_two_level_velocity_solver,
                    E->control.ala_two_level_velocity_iterations,
                    E->control.ala_two_level_velocity_eigenvalue_min,
                    E->control.ala_two_level_velocity_eigenvalue_max,
                    global_invalid_coarse);
        }
    }
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=0; j<neq; j++)
            F[m][j] = FF[m][j];

    /* FF contains the current -C^T*P forcing.  Pressure increments below
     * apply the complete strict-ALA transpose explicitly. */
    initial_vel_residual(E, V, P, F, imp);

    assemble_div_rho_u(E, V, r, lev);
    residual = incompressibility_residual(E, V, r);
    initial_rnorm = sqrt(global_pdot(E, r, r, lev));
    relative_residual = (initial_rnorm > 0.0) ? 1.0 : 0.0;
    strict_ala_continuity_metrics(E, V, r, div_u, lev,
                                  &initial_mass_norm, &cancellation_l2);
    mass_norm = initial_mass_norm;
    mass_relative_residual = (initial_mass_norm > 0.0) ? 1.0 : 0.0;

    sq_vdotv = sqrt(E->monitor.vdotv);
    dpressure = 1.0;
    dvelocity = 1.0;
    if(E->control.print_convergence && E->parallel.me == 0) {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
        fprintf(E->fp,
                "ALA PCG continuity residuals: cancellation=%e "
                "mass_relative=%e algebraic_relative=%e\n",
                cancellation_l2, mass_relative_residual,
                relative_residual);
    }
    strict_ala_depth_diagnostics(E, r, div_u, lev, count);
    strict_ala_coarse_residual_diagnostics(E, r, lev, count);

    /* A fixed inner target keeps the approximate Schur operator as
     * stationary as the RHS-relative multigrid stopping rule permits. */
    inner_relative_accuracy = E->control.ala_inner_accuracy_max;
    rho_old = 0.0;
    restart_search = 1;
    hybrid_consecutive_count = 0;
    hybrid_converged = 0;
    nonpositive_curvature_count = 0;
    min_curvature = 1.0e300;

    while(count < *steps_max &&
          (cancellation_l2 >= E->control.tole_comp ||
           (E->control.ala_hybrid_convergence && !hybrid_converged))) {

        apply_ala_pressure_preconditioner(
            E,r,z,preconditioner_work,lev,count,&patch_cache);

        rho = global_pdot(E, r, z, lev);
        if(!isfinite(rho) || rho <= 1.0e-300) {
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA PCG breakdown: invalid preconditioned residual "
                        "<r,z>=%e\n", rho);
                fprintf(stderr,
                        "ALA PCG breakdown: invalid preconditioned residual "
                        "<r,z>=%e\n", rho);
            }
            parallel_process_termination();
        }

        if(restart_search) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p[m][j] = z[m][j];
            restart_search = 0;
        }
        else {
            beta = rho / rho_old;
            if(!isfinite(beta)) {
                if(E->parallel.me == 0) {
                    fprintf(E->fp, "ALA PCG breakdown: non-finite beta\n");
                    fprintf(stderr, "ALA PCG breakdown: non-finite beta\n");
                }
                parallel_process_termination();
            }
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p[m][j] = z[m][j] + beta * p[m][j];
        }

        assemble_grad_rho_p(E, p, F, lev);
        inner_accuracy = strict_ala_inner_accuracy(
            E, F, lev, inner_relative_accuracy);
        valid = solve_del2_u(E, E->u1, F, inner_accuracy, lev);
        if(!valid) {
            if(E->parallel.me == 0) {
                fputs("ALA PCG inner velocity solve failed\n", stderr);
                fputs("ALA PCG inner velocity solve failed\n", E->fp);
            }
            parallel_process_termination();
        }
        strip_bcs_from_residual(E, E->u1, lev);
        assemble_div_rho_u(E, E->u1, q, lev);

        curvature = global_pdot(E, p, q, lev);
        if(!isfinite(curvature) || curvature <= 1.0e-300) {
            nonpositive_curvature_count++;
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA PCG breakdown: nonpositive Schur curvature=%e "
                        "at iteration %d\n", curvature, count);
                fprintf(stderr,
                        "ALA PCG breakdown: nonpositive Schur curvature=%e "
                        "at iteration %d\n", curvature, count);
            }
            parallel_process_termination();
        }
        min_curvature = min(min_curvature, curvature);
        alpha = rho / curvature;
        if(!isfinite(alpha)) {
            if(E->parallel.me == 0) {
                fprintf(E->fp, "ALA PCG breakdown: non-finite alpha\n");
                fprintf(stderr, "ALA PCG breakdown: non-finite alpha\n");
            }
            parallel_process_termination();
        }

        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            for(j=1; j<=npno; j++) {
                P[m][j] += alpha * p[m][j];
                r[m][j] -= alpha * q[m][j];
            }
            for(j=0; j<neq; j++)
                V[m][j] -= alpha * E->u1[m][j];
        }

        recursive_rnorm = sqrt(global_pdot(E, r, r, lev));
        recursive_relative_residual = (initial_rnorm > 0.0)
            ? recursive_rnorm / initial_rnorm : 0.0;

        assemble_div_rho_u(E, V, explicit_r, lev);
        residual = incompressibility_residual(E, V, explicit_r);
        rnorm = sqrt(global_pdot(E, explicit_r, explicit_r, lev));
        relative_residual = (initial_rnorm > 0.0)
            ? rnorm / initial_rnorm : 0.0;
        strict_ala_continuity_metrics(E, V, explicit_r, div_u, lev,
                                      &mass_norm, &cancellation_l2);
        mass_relative_residual = (initial_mass_norm > 0.0)
            ? mass_norm / initial_mass_norm : 0.0;
        drift_ratio = (relative_residual > 0.0)
            ? recursive_relative_residual / relative_residual : 1.0;
        residual = cancellation_l2;

        dpressure = fabs(alpha) * sqrt(
            global_pdot(E, p, p, lev) /
            (1.0e-32 + global_pdot(E, P, P, lev)));
        dvelocity = fabs(alpha) * sqrt(
            global_vdot(E, E->u1, E->u1, lev) /
            (1.0e-32 + E->monitor.vdotv));

        count++;
        if(E->control.ala_hybrid_convergence) {
            if(E->monitor.incompressibility <
                   E->control.ala_div_v_tolerance &&
               dvelocity < E->control.ala_update_tolerance &&
               dpressure < E->control.ala_update_tolerance)
                hybrid_consecutive_count++;
            else
                hybrid_consecutive_count = 0;
            if(hybrid_consecutive_count >= E->control.ala_consecutive_steps)
                hybrid_converged = 1;
        }

        sq_vdotv = sqrt(E->monitor.vdotv);
        if(E->control.print_convergence && E->parallel.me == 0) {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
            fprintf(E->fp,
                    "ALA PCG continuity residuals: cancellation=%e "
                    "mass_relative=%e algebraic_relative=%e "
                    "recursive_algebraic=%e drift=%e inner_rel=%e "
                    "curvature=%e\n",
                    cancellation_l2, mass_relative_residual,
                    relative_residual, recursive_relative_residual,
                    drift_ratio, inner_relative_accuracy, curvature);
            if(E->control.ala_hybrid_convergence)
                fprintf(E->fp,
                        "ALA hybrid convergence streak = %d/%d "
                        "limits: div/v<%e updates<%e\n",
                        hybrid_consecutive_count,
                        E->control.ala_consecutive_steps,
                        E->control.ala_div_v_tolerance,
                        E->control.ala_update_tolerance);
        }
        strict_ala_depth_diagnostics(E, explicit_r, div_u, lev, count);
        strict_ala_coarse_residual_diagnostics(E, explicit_r, lev, count);

        rho_old = rho;

        /* Explicit replacement protects a long inexact PCG run from residual
         * drift.  Replacement also restarts the conjugate direction. */
        if((count % E->control.ala_pcg_restart_interval) == 0 ||
           drift_ratio > 10.0 || drift_ratio < 0.1) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    r[m][j] = explicit_r[m][j];
            restart_search = 1;
            if(E->parallel.me == 0)
                fprintf(E->fp,
                        "ALA PCG restarted from explicit residual at "
                        "iteration %d\n", count);
        }
    }

    if(E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA PCG operator audit: minimum_curvature=%e "
                "nonpositive_curvature=%d iterations=%d "
                "Schur_applications=%d pressure_velocity_solves=%d\n",
                min_curvature, nonpositive_curvature_count, count, count,
                count);
        fprintf(stderr,
                "ALA PCG operator audit: minimum_curvature=%e "
                "nonpositive_curvature=%d iterations=%d "
                "Schur_applications=%d pressure_velocity_solves=%d\n",
                min_curvature, nonpositive_curvature_count, count, count,
                count);
    }

    if(E->parallel.me == 0 &&
       cancellation_l2 < E->control.tole_comp &&
       (!E->control.ala_hybrid_convergence || hybrid_converged)) {
        fprintf(E->fp,
                "Strict ALA PCG converged by physical continuity: "
                "cancellation=%e tolerance=%e mass_relative=%e "
                "algebraic_relative=%e iterations=%d\n",
                cancellation_l2, E->control.tole_comp,
                mass_relative_residual, relative_residual, count);
        fflush(E->fp);
    }

    if(cancellation_l2 >= E->control.tole_comp ||
       (E->control.ala_hybrid_convergence && !hybrid_converged)) {
        if(E->parallel.me == 0) {
            fprintf(E->fp,
                    "Strict ALA PCG failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d "
                    "iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fprintf(stderr,
                    "Strict ALA PCG failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d "
                    "iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fflush(E->fp);
        }
        parallel_process_termination();
    }

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free((void *)F[m]);
        free((void *)r[m]);
        free((void *)z[m]);
        free((void *)p[m]);
        free((void *)q[m]);
        free((void *)explicit_r[m]);
        free((void *)div_u[m]);
        free((void *)preconditioner_work[m]);
    }
    if(E->control.ala_shallow_patch_preconditioner)
        free_ala_shallow_patch_cache(E,&patch_cache);

    *steps_max = count;
    return(residual);
}


/* Solve compressible Stokes flow using
 * bi-conjugate gradient stablized (BiCG-stab) iterations
 */

static float solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **FF,
                                    double imp, int *steps_max)
{
    void assemble_div_rho_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void assemble_grad_rho_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();

    double global_vdot(), global_pdot();
    double CPU_time0();

    int gnpno, gneq;
    int npno, neq;
    int m, j, count, lev;
    int valid;
    int restart_search;
    int hybrid_consecutive_count, hybrid_converged;

    double alpha, beta, omega,sq_vdotv;
    double r0dotrt, r1dotrt;
    double residual, dpressure, dvelocity;
    double initial_rnorm, relative_residual, rnorm;
    double initial_mass_norm, mass_norm, mass_relative_residual;
    double cancellation_l2;
    double recursive_rnorm, recursive_relative_residual, drift_ratio;
    double denominator, numerator, inner_accuracy;
    double inner_relative_accuracy, previous_inner_relative_accuracy;
    double symmetry_left, symmetry_right, symmetry_defect;
    double max_symmetry_defect, curvature_p, curvature_s;
    int symmetry_sample_count, nonpositive_curvature_count;

    double *F[NCS];
    double *r1[NCS], *r2[NCS], *pt[NCS], *p1[NCS], *p2[NCS];
    double *rt[NCS], *v0[NCS], *s0[NCS], *st[NCS], *t0[NCS];
    double *u0[NCS], *div_u[NCS];
    double *shuffle[NCS];

    double time0, v_res;

    gnpno = E->mesh.npno;
    gneq = E->mesh.neq;
    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[m] = (double *)malloc(neq*sizeof(double));
        r1[m] = (double *)malloc((npno+1)*sizeof(double));
        r2[m] = (double *)malloc((npno+1)*sizeof(double));
        pt[m] = (double *)malloc((npno+1)*sizeof(double));
        p1[m] = (double *)malloc((npno+1)*sizeof(double));
        p2[m] = (double *)malloc((npno+1)*sizeof(double));
        rt[m] = (double *)malloc((npno+1)*sizeof(double));
        v0[m] = (double *)malloc((npno+1)*sizeof(double));
        s0[m] = (double *)malloc((npno+1)*sizeof(double));
        st[m] = (double *)malloc((npno+1)*sizeof(double));
        t0[m] = (double *)malloc((npno+1)*sizeof(double));
        div_u[m] = (double *)malloc((npno+1)*sizeof(double));

        u0[m] = (double *)malloc(neq*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[m][j] = FF[m][j];


    /* FF already contains the current -C^T*P force contribution.  The
     * pressure increment operator below must therefore apply both the
     * standard gradient and the strict-ALA C^T response. */
    /* calculate the initial velocity residual */
    v_res = initial_vel_residual(E, V, P, F, imp);


    /* initial residual r1 = div(rho_ref*V) */
    assemble_div_rho_u(E, V, r1, lev);
    residual = incompressibility_residual(E, V, r1);
    initial_rnorm = sqrt(global_pdot(E, r1, r1, lev));
    relative_residual = (initial_rnorm > 0.0) ? 1.0 : 0.0;
    if(E->control.ala_pressure_buoyancy) {
        strict_ala_continuity_metrics(E, V, r1, div_u, lev,
                                      &initial_mass_norm, &cancellation_l2);
        mass_norm = initial_mass_norm;
        mass_relative_residual = (initial_mass_norm > 0.0) ? 1.0 : 0.0;
    }
    else {
        initial_mass_norm = mass_norm = initial_rnorm;
        mass_relative_residual = relative_residual;
        cancellation_l2 = relative_residual;
    }


    /* initial conjugate residual rt = r1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=1; j<=npno; j++)
            rt[m][j] = r1[m][j];


    sq_vdotv = sqrt(E->monitor.vdotv);

    /* pressure and velocity corrections */
    dpressure = 1.0;
    dvelocity = 1.0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0, sq_vdotv,
                                   dvelocity, dpressure);
        if(E->control.ala_pressure_buoyancy)
            fprintf(E->fp,
                    "ALA continuity residuals: cancellation=%e "
                    "mass_relative=%e algebraic_relative=%e\n",
                    cancellation_l2, mass_relative_residual,
                    relative_residual);
    }
    if(E->control.ala_pressure_buoyancy)
        strict_ala_depth_diagnostics(E, r1, div_u, lev, count);


    valid = 1;
    r0dotrt = alpha = omega = 0;
    restart_search = 0;
    previous_inner_relative_accuracy = -1.0;
    max_symmetry_defect = 0.0;
    symmetry_sample_count = 0;
    nonpositive_curvature_count = 0;
    hybrid_consecutive_count = 0;
    hybrid_converged = 0;

    while( (count < *steps_max) &&
           ((E->control.ala_pressure_buoyancy &&
             (cancellation_l2 >= E->control.tole_comp ||
              (E->control.ala_hybrid_convergence &&
               !hybrid_converged))) ||
            (!E->control.ala_pressure_buoyancy &&
             E->monitor.incompressibility >= E->control.tole_comp)) &&
           (E->control.ala_pressure_buoyancy ||
            ((dpressure >= imp) && (dvelocity >= imp))) )  {


        inner_relative_accuracy = strict_ala_inner_relative_accuracy(
            E, mass_relative_residual);
        if(previous_inner_relative_accuracy > 0.0 &&
           fabs(inner_relative_accuracy - previous_inner_relative_accuracy)
             > 1.0e-12 * previous_inner_relative_accuracy) {
            /* A changed inner tolerance changes the approximate K inverse.
             * Re-anchor BiCGStab before using the new inexact Schur operator. */
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++) {
                    r1[m][j] = t0[m][j];
                    rt[m][j] = t0[m][j];
                }
            restart_search = 1;
            r0dotrt = alpha = omega = 0.0;
            if(E->parallel.me == 0) {
                fprintf(E->fp,
                        "ALA inner relative accuracy changed %.3e -> %.3e; "
                        "restarting BiCGStab at iteration %d\n",
                        previous_inner_relative_accuracy,
                        inner_relative_accuracy, count);
                fprintf(stderr,
                        "ALA inner relative accuracy changed %.3e -> %.3e; "
                        "restarting BiCGStab at iteration %d\n",
                        previous_inner_relative_accuracy,
                        inner_relative_accuracy, count);
            }
        }
        previous_inner_relative_accuracy = inner_relative_accuracy;


        /* r1dotrt = <r1, rt> */
        r1dotrt = global_pdot(E, r1, rt, lev);
        if(!isfinite(r1dotrt) || fabs(r1dotrt) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid <r,rt>\n");
            fprintf(stderr, "BiCGStab breakdown: invalid <r,rt>\n");
            parallel_process_termination();
        }


        /* update search direction */
        if(count == 0 || restart_search) {
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j];
            restart_search = 0;
        }
        else {
            /* p2 = r1 + <r1,rt>/<r0,rt> * alpha/omega * (p1 - omega*v0) */
            if(!isfinite(r0dotrt) || fabs(r0dotrt) <= 1.0e-300 ||
               !isfinite(omega) || fabs(omega) <= 1.0e-300) {
                fprintf(E->fp,
                        "BiCGStab breakdown: invalid beta denominator\n");
                fprintf(stderr,
                        "BiCGStab breakdown: invalid beta denominator\n");
                parallel_process_termination();
            }
            beta = (r1dotrt / r0dotrt) * (alpha / omega);
            if(!isfinite(beta)) {
                fprintf(E->fp, "BiCGStab breakdown: non-finite beta\n");
                fprintf(stderr, "BiCGStab breakdown: non-finite beta\n");
                parallel_process_termination();
            }
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++)
                    p2[m][j] = r1[m][j] + beta
                        * (p1[m][j] - omega * v0[m][j]);
        }


        /* pressure Schur-diagonal preconditioner */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                pt[m][j] = E->BPI[lev][m][j] * p2[m][j];


        /* solve K*u0 = (D+C)^T*pt for strict ALA */
        if(E->control.ala_pressure_buoyancy)
            assemble_grad_rho_p(E, pt, F, lev);
        else
            assemble_grad_p(E, pt, F, lev);
        inner_accuracy = E->control.ala_pressure_buoyancy
            ? strict_ala_inner_accuracy(E, F, lev,
                                        inner_relative_accuracy) : imp * v_res;
        valid = solve_del2_u(E, u0, F, inner_accuracy, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, u0, lev);


        /* v0 = div(rho_ref*u0) */
        assemble_div_rho_u(E, u0, v0, lev);


        /* alpha = r1dotrt / <rt, v0> */
        denominator = global_pdot(E, rt, v0, lev);
        if(!isfinite(denominator) || fabs(denominator) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid alpha denominator\n");
            fprintf(stderr, "BiCGStab breakdown: invalid alpha denominator\n");
            parallel_process_termination();
        }
        alpha = r1dotrt / denominator;
        if(!isfinite(alpha)) {
            fprintf(E->fp, "BiCGStab breakdown: non-finite alpha\n");
            fprintf(stderr, "BiCGStab breakdown: non-finite alpha\n");
            parallel_process_termination();
        }


        /* s0 = r1 - alpha * v0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = r1[m][j] - alpha * v0[m][j];


        /* pressure Schur-diagonal preconditioner */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                st[m][j] = E->BPI[lev][m][j] * s0[m][j];


        /* solve K*u1 = (D+C)^T*st for strict ALA */
        if(E->control.ala_pressure_buoyancy)
            assemble_grad_rho_p(E, st, F, lev);
        else
            assemble_grad_p(E, st, F, lev);
        inner_accuracy = E->control.ala_pressure_buoyancy
            ? strict_ala_inner_accuracy(E, F, lev,
                                        inner_relative_accuracy) : imp * v_res;
        valid = solve_del2_u(E, E->u1, F, inner_accuracy, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 2\n", stderr);
            fputs("Warning: solver not converging! 2\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* t0 = div(rho_ref * u1) */
        assemble_div_rho_u(E, E->u1, t0, lev);

        if(E->control.ala_pressure_buoyancy &&
           E->control.ala_schur_symmetry_check) {
            /* v0=S*pt and t0=S*st are already available.  This checks the
             * actual inexact Schur applications without another K solve. */
            symmetry_left = global_pdot(E, pt, t0, lev);
            symmetry_right = global_pdot(E, st, v0, lev);
            symmetry_defect = fabs(symmetry_left - symmetry_right)
                / (fabs(symmetry_left) + fabs(symmetry_right) + 1.0e-300);
            curvature_p = global_pdot(E, pt, v0, lev);
            curvature_s = global_pdot(E, st, t0, lev);
            if(!isfinite(symmetry_defect) || !isfinite(curvature_p) ||
               !isfinite(curvature_s)) {
                symmetry_defect = 1.0;
                nonpositive_curvature_count++;
            }
            else if(curvature_p <= 0.0 || curvature_s <= 0.0)
                nonpositive_curvature_count++;
            max_symmetry_defect = max(max_symmetry_defect, symmetry_defect);
            symmetry_sample_count++;
        }


        /* omega = <t0, s0> / <t0, t0> */
        numerator = global_pdot(E, t0, s0, lev);
        denominator = global_pdot(E, t0, t0, lev);
        if(!isfinite(numerator) || !isfinite(denominator) ||
           fabs(denominator) <= 1.0e-300) {
            fprintf(E->fp, "BiCGStab breakdown: invalid omega denominator\n");
            fprintf(stderr, "BiCGStab breakdown: invalid omega denominator\n");
            parallel_process_termination();
        }
        omega = numerator / denominator;
        if(!isfinite(omega)) {
            fprintf(E->fp, "BiCGStab breakdown: non-finite omega\n");
            fprintf(stderr, "BiCGStab breakdown: non-finite omega\n");
            parallel_process_termination();
        }


        /* r2 = s0 - omega * t0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                r2[m][j] = s0[m][j] - omega * t0[m][j];

        recursive_rnorm = sqrt(global_pdot(E, r2, r2, lev));
        recursive_relative_residual = (initial_rnorm > 0.0)
            ? recursive_rnorm / initial_rnorm : 0.0;


        /* P = P + alpha * pt + omega * st */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                s0[m][j] = alpha * pt[m][j] + omega * st[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=1; j<=npno; j++)
                P[m][j] += s0[m][j];


        /* V = V - alpha * u0 - omega * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                F[m][j] = alpha * u0[m][j] + omega * E->u1[m][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[m][j] -= F[m][j];


        /* compute velocity and incompressibility residual */
        assemble_div_rho_u(E, V, t0, lev);
        residual = incompressibility_residual(E, V, t0);
        rnorm = sqrt(global_pdot(E, t0, t0, lev));
        relative_residual = (initial_rnorm > 0.0)
            ? rnorm / initial_rnorm : 0.0;
        if(E->control.ala_pressure_buoyancy) {
            strict_ala_continuity_metrics(E, V, t0, div_u, lev,
                                          &mass_norm, &cancellation_l2);
            mass_relative_residual = (initial_mass_norm > 0.0)
                ? mass_norm / initial_mass_norm : 0.0;
        }
        drift_ratio = (relative_residual > 0.0)
            ? recursive_relative_residual / relative_residual : 1.0;
        if(E->control.ala_pressure_buoyancy)
            residual = cancellation_l2;

        /* compute velocity and pressure corrections */
        dpressure = sqrt( global_pdot(E, s0, s0, lev)
                          / (1.0e-32 + global_pdot(E, P, P, lev)) );
        dvelocity = sqrt( global_vdot(E, F, F, lev)
                          / (1.0e-32 + E->monitor.vdotv) );


        count++;

        if(E->control.ala_pressure_buoyancy &&
           E->control.ala_hybrid_convergence) {
            if(E->monitor.incompressibility <
                   E->control.ala_div_v_tolerance &&
               dvelocity < E->control.ala_update_tolerance &&
               dpressure < E->control.ala_update_tolerance)
                hybrid_consecutive_count++;
            else
                hybrid_consecutive_count = 0;
            if(hybrid_consecutive_count >= E->control.ala_consecutive_steps)
                hybrid_converged = 1;
        }

	sq_vdotv = sqrt(E->monitor.vdotv);

        if(E->control.print_convergence && E->parallel.me==0) {
            print_convergence_progress(E, count, time0, sq_vdotv,
                                       dvelocity, dpressure);
            if(E->control.ala_pressure_buoyancy)
                fprintf(E->fp,
                        "ALA continuity residuals: cancellation=%e "
                        "mass_relative=%e algebraic_relative=%e "
                        "recursive_algebraic=%e drift=%e inner_rel=%e\n",
                        cancellation_l2, mass_relative_residual,
                        relative_residual, recursive_relative_residual,
                        drift_ratio, inner_relative_accuracy);
            if(E->control.ala_pressure_buoyancy &&
               E->control.ala_hybrid_convergence)
                fprintf(E->fp,
                        "ALA hybrid convergence streak = %d/%d "
                        "limits: div/v<%e updates<%e\n",
                        hybrid_consecutive_count,
                        E->control.ala_consecutive_steps,
                        E->control.ala_div_v_tolerance,
                        E->control.ala_update_tolerance);
            if(E->control.ala_pressure_buoyancy &&
               E->control.ala_schur_symmetry_check) {
                fprintf(E->fp,
                        "ALA Schur symmetry defect = %e curvature=(%e,%e)\n",
                        symmetry_defect, curvature_p, curvature_s);
                fprintf(stderr,
                        "ALA Schur symmetry defect = %e curvature=(%e,%e)\n",
                        symmetry_defect, curvature_p, curvature_s);
            }
        }
        if(E->control.ala_pressure_buoyancy)
            strict_ala_depth_diagnostics(E, t0, div_u, lev, count);

        /* Re-anchor the recurrence to the explicitly assembled B*u residual.
         * This limits finite-precision drift in long strict-ALA solves. */
        if(E->control.ala_pressure_buoyancy &&
           ((count % 20) == 0 || drift_ratio > 10.0 || drift_ratio < 0.1)) {
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=1; j<=npno; j++) {
                    r2[m][j] = t0[m][j];
                    rt[m][j] = t0[m][j];
                }
            restart_search = 1;
            if(E->parallel.me == 0)
                fprintf(E->fp,
                        "ALA BiCGStab restarted from explicit residual at "
                        "iteration %d\n", count);
        }


        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[m] = p1[m];
            p1[m] = p2[m];
            p2[m] = shuffle[m];

            shuffle[m] = r1[m];
            r1[m] = r2[m];
            r2[m] = shuffle[m];
        }

        /* shift <r0, rt> = <r1, rt> */
        r0dotrt = r1dotrt;

    } /* end loop for BiCGStab */

    if(E->control.ala_pressure_buoyancy &&
       E->control.ala_schur_symmetry_check && E->parallel.me == 0) {
        fprintf(E->fp,
                "ALA Schur CG suitability: max_symmetry_defect=%e "
                "tolerance=%e samples=%d nonpositive_curvature=%d result=%s\n",
                max_symmetry_defect,
                E->control.ala_schur_symmetry_tolerance,
                symmetry_sample_count,
                nonpositive_curvature_count,
                (symmetry_sample_count > 0 && max_symmetry_defect <=
                     E->control.ala_schur_symmetry_tolerance &&
                 nonpositive_curvature_count == 0) ? "PASS" : "FAIL");
        fprintf(stderr,
                "ALA Schur CG suitability: max_symmetry_defect=%e "
                "tolerance=%e samples=%d nonpositive_curvature=%d result=%s\n",
                max_symmetry_defect,
                E->control.ala_schur_symmetry_tolerance,
                symmetry_sample_count,
                nonpositive_curvature_count,
                (symmetry_sample_count > 0 && max_symmetry_defect <=
                     E->control.ala_schur_symmetry_tolerance &&
                 nonpositive_curvature_count == 0) ? "PASS" : "FAIL");
    }

    if(E->control.ala_pressure_buoyancy && E->parallel.me == 0 &&
       cancellation_l2 < E->control.tole_comp &&
       (!E->control.ala_hybrid_convergence || hybrid_converged)) {
        fprintf(E->fp,
                "Strict ALA BiCGStab converged by physical continuity: "
                "cancellation=%e tolerance=%e mass_relative=%e "
                "algebraic_relative=%e iterations=%d\n",
                cancellation_l2, E->control.tole_comp,
                mass_relative_residual, relative_residual, count);
        fflush(E->fp);
    }

    if(E->control.ala_pressure_buoyancy &&
       (cancellation_l2 >= E->control.tole_comp ||
        (E->control.ala_hybrid_convergence && !hybrid_converged))) {
        if(E->parallel.me == 0) {
            fprintf(stderr,
                    "Strict ALA BiCGStab failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fprintf(E->fp,
                    "Strict ALA BiCGStab failed physical continuity: "
                    "cancellation=%e tolerance=%e mass_relative=%e "
                    "algebraic_relative=%e hybrid_streak=%d/%d iterations=%d\n",
                    cancellation_l2, E->control.tole_comp,
                    mass_relative_residual, relative_residual,
                    hybrid_consecutive_count,
                    E->control.ala_consecutive_steps, count);
            fflush(E->fp);
        }
        parallel_process_termination();
    }


    for(m=1; m<=E->sphere.caps_per_proc; m++) {
    	free((void *) F[m]);
        free((void *) r1[m]);
        free((void *) r2[m]);
        free((void *) pt[m]);
        free((void *) p1[m]);
        free((void *) p2[m]);
        free((void *) rt[m]);
        free((void *) v0[m]);
        free((void *) s0[m]);
        free((void *) st[m]);
        free((void *) t0[m]);
        free((void *) div_u[m]);

        free((void *) u0[m]);
    }

    *steps_max=count;

    return(residual);

}


/* Solve compressible Stokes flow using
 * conjugate gradient (CG) iterations with an outer iteration
 */

static float solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max)
{
    int m, i;
    int cycles, num_of_loop;
    double residual;
    double relative_err_v, relative_err_p;
    double *old_v[NCS], *old_p[NCS],*diff_v[NCS],*diff_p[NCS];

    const int npno = E->lmesh.npno;
    const int neq = E->lmesh.neq;
    const int lev = E->mesh.levmax;

    double global_vdot(),global_pdot();
    void assemble_forces();

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	old_v[m] = (double *)malloc(neq*sizeof(double));
    	diff_v[m] = (double *)malloc(neq*sizeof(double));
    	old_p[m] = (double *)malloc((npno+1)*sizeof(double));
    	diff_p[m] = (double *)malloc((npno+1)*sizeof(double));
    }

    cycles = E->control.p_iterations;

    residual = 1.0;
    relative_err_v = 1.0;
    relative_err_p = 1.0;
    num_of_loop = 0;

    while((relative_err_v >= E->control.relative_err_accuracy ||
           relative_err_p >= E->control.relative_err_accuracy) &&
          num_of_loop <= E->control.compress_iter_maxstep) {

        for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++) old_v[m][i] = V[m][i];
            for(i=1;i<=npno;i++) old_p[m][i] = P[m][i];
        }

        if(E->control.ala_pressure_buoyancy)
            assemble_forces(E,0);

        residual = solve_Ahat_p_fhat_CG(E,V,P,F,E->control.accuracy,&cycles);

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=0;i<neq;i++) diff_v[m][i] = V[m][i] - old_v[m][i];

        relative_err_v = sqrt( global_vdot(E,diff_v,diff_v,lev) /
                               (1.0e-32 + global_vdot(E,V,V,lev)) );

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=npno;i++) diff_p[m][i] = P[m][i] - old_p[m][i];

        relative_err_p = sqrt( global_pdot(E,diff_p,diff_p,lev) /
                               (1.0e-32 + global_pdot(E,P,P,lev)) );

        num_of_loop++;

        if(E->parallel.me == 0) {
            fprintf(stderr, "Relative error err_v / v = %e and err_p / p = %e after %d loops\n\n", relative_err_v, relative_err_p, num_of_loop);
            fprintf(E->fp, "Relative error err_v / v = %e and err_p / p = %e after %d loops\n\n", relative_err_v, relative_err_p, num_of_loop);
        }

    } /* end of while */

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	free((void *) old_v[m]);
    	free((void *) old_p[m]);
	free((void *) diff_v[m]);
	free((void *) diff_p[m]);
    }

    return(residual);
}


static double initial_vel_residual(struct All_variables *E,
                                   double **V, double **P, double **F,
                                   double imp)
{
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    double global_vdot();

    int neq = E->lmesh.neq;
    int gneq = E->mesh.neq;
    int lev = E->mesh.levmax;
    int i, m, valid;
    double v_res;

    v_res = sqrt(global_vdot(E, F, F, lev) / gneq);

    if (E->parallel.me==0) {
        fprintf(E->fp, "initial residue of momentum equation F %.9e %d\n",
                v_res, gneq);
        fprintf(stderr, "initial residue of momentum equation F %.9e %d\n",
                v_res, gneq);
    }


    /* F = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    assemble_del2_u(E, V, E->u1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[m][i] = F[m][i] - E->u1[m][i];

    strip_bcs_from_residual(E, F, lev);


    /* solve K*u1 = F for u1 */
    valid = solve_del2_u(E, E->u1, F, imp*v_res, lev);
    if(!valid && (E->parallel.me==0)) {
        fputs("Warning: solver not converging! 0\n", stderr);
        fputs("Warning: solver not converging! 0\n", E->fp);
    }
    strip_bcs_from_residual(E, E->u1, lev);


    /* V = V + u1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            V[m][i] += E->u1[m][i];

    return(v_res);
}



static double incompressibility_residual(struct All_variables *E,
                                         double **V, double **r)
{
    double global_pdot();
    double global_vdot();

    int gnpno = E->mesh.npno;
    int gneq = E->mesh.neq;
    int lev = E->mesh.levmax;
    double tmp1, tmp2;

    /* incompressiblity residual = norm(r) / norm(V) */

    tmp1 = global_vdot(E, V, V, lev);
    tmp2 = global_pdot(E, r, r, lev);
    E->monitor.incompressibility = sqrt((gneq / gnpno)
                                        *( (1.0e-32 + tmp2)
                                              / (1.0e-32 + tmp1) ));

    E->monitor.vdotv = tmp1;

    return(sqrt(tmp2/gnpno));;
}
