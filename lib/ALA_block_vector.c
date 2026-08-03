/* Distributed vector algebra for a future monolithic strict-ALA Krylov path. */

#include <math.h>
#include <stdlib.h>

#include "global_defs.h"
#include "ala_block_vector.h"

void myerror(struct All_variables *,char *);


static void ala_block_vector_require_level(
    struct All_variables *E, const struct ala_block_vector *vector, int level)
{
    if(vector==NULL || vector->velocity==NULL || vector->pressure==NULL)
        myerror(E,"Invalid strict-ALA block vector");
    if(vector->level!=level)
        myerror(E,"Mismatched level in strict-ALA block vector operation");
}


static double **ala_block_allocate_field(struct All_variables *E, int entries)
{
    int m;
    double **field;

    field=(double **)calloc(NCS,sizeof(double *));
    if(field==NULL)
        myerror(E,"Unable to allocate strict-ALA block-vector pointers");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        field[m]=(double *)calloc(entries,sizeof(double));
        if(field[m]==NULL)
            myerror(E,"Unable to allocate strict-ALA block-vector field");
    }
    return field;
}


struct ala_block_vector *ala_block_vector_create(struct All_variables *E,
                                                 int level)
{
    struct ala_block_vector *vector;
    int neq,npno;

    if(level<E->mesh.levmin || level>E->mesh.levmax)
        myerror(E,"Invalid level for strict-ALA block vector");
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    vector=(struct ala_block_vector *)calloc(1,sizeof(*vector));
    if(vector==NULL)
        myerror(E,"Unable to allocate strict-ALA block vector");
    vector->velocity=ala_block_allocate_field(E,neq+1);
    vector->pressure=ala_block_allocate_field(E,npno+1);
    vector->level=level;
    return vector;
}


void ala_block_vector_destroy(struct All_variables *E,
                              struct ala_block_vector *vector)
{
    int m;

    if(vector==NULL)
        return;
    if(vector->velocity!=NULL) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            free(vector->velocity[m]);
        free(vector->velocity);
    }
    if(vector->pressure!=NULL) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            free(vector->pressure[m]);
        free(vector->pressure);
    }
    free(vector);
}


void ala_block_vector_zero(struct All_variables *E,
                           struct ala_block_vector *vector)
{
    int m,i,neq,npno,level;

    if(vector==NULL)
        myerror(E,"Invalid strict-ALA block vector");
    level=vector->level;
    ala_block_vector_require_level(E,vector,level);
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<=neq;i++)
            vector->velocity[m][i]=0.0;
        for(i=0;i<=npno;i++)
            vector->pressure[m][i]=0.0;
    }
}


void ala_block_vector_copy(struct All_variables *E,
                           const struct ala_block_vector *source,
                           struct ala_block_vector *destination)
{
    int m,i,neq,npno,level;

    if(source==NULL || destination==NULL)
        myerror(E,"Invalid strict-ALA block vector copy");
    level=source->level;
    ala_block_vector_require_level(E,source,level);
    ala_block_vector_require_level(E,destination,level);
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            destination->velocity[m][i]=source->velocity[m][i];
        destination->velocity[m][neq]=0.0;
        destination->pressure[m][0]=0.0;
        for(i=1;i<=npno;i++)
            destination->pressure[m][i]=source->pressure[m][i];
    }
}


void ala_block_vector_scale(struct All_variables *E, double scale,
                            struct ala_block_vector *vector)
{
    int m,i,neq,npno,level;

    if(!isfinite(scale))
        myerror(E,"Non-finite strict-ALA block-vector scale");
    if(vector==NULL)
        myerror(E,"Invalid strict-ALA block vector scale");
    level=vector->level;
    ala_block_vector_require_level(E,vector,level);
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            vector->velocity[m][i]*=scale;
        vector->velocity[m][neq]=0.0;
        vector->pressure[m][0]=0.0;
        for(i=1;i<=npno;i++)
            vector->pressure[m][i]*=scale;
    }
}


void ala_block_vector_axpy(struct All_variables *E, double scale,
                           const struct ala_block_vector *source,
                           struct ala_block_vector *destination)
{
    int m,i,neq,npno,level;

    if(!isfinite(scale))
        myerror(E,"Non-finite strict-ALA block-vector axpy scale");
    if(source==NULL || destination==NULL)
        myerror(E,"Invalid strict-ALA block vector axpy");
    level=source->level;
    ala_block_vector_require_level(E,source,level);
    ala_block_vector_require_level(E,destination,level);
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            destination->velocity[m][i] += scale*source->velocity[m][i];
        destination->velocity[m][neq]=0.0;
        destination->pressure[m][0]=0.0;
        for(i=1;i<=npno;i++)
            destination->pressure[m][i] += scale*source->pressure[m][i];
    }
}


static void ala_block_vector_component_dot(
    struct All_variables *E, const struct ala_block_vector *left,
    const struct ala_block_vector *right, double component_dot[2])
{
    int m,i,level,neq,npno;
    double local[2],global[2],duplicate;

    if(left==NULL || right==NULL)
        myerror(E,"Invalid strict-ALA block-vector dot product");
    level=left->level;
    ala_block_vector_require_level(E,left,level);
    ala_block_vector_require_level(E,right,level);
    neq=E->lmesh.NEQ[level];
    npno=E->lmesh.NPNO[level];
    local[0]=0.0;
    local[1]=0.0;
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        duplicate=0.0;
        for(i=0;i<neq;i++)
            local[0] += left->velocity[m][i]*right->velocity[m][i];
        for(i=1;i<=E->parallel.Skip_neq[level][m];i++)
            duplicate += left->velocity[m][E->parallel.Skip_id[level][m][i]]
                       * right->velocity[m][E->parallel.Skip_id[level][m][i]];
        local[0] -= duplicate;
        for(i=1;i<=npno;i++) {
            if(E->ECO[level][m][i].area<=0.0)
                myerror(E,"Nonpositive pressure volume in strict-ALA "
                          "block residual metric");
            local[1] += left->pressure[m][i]*right->pressure[m][i]
                       /E->ECO[level][m][i].area;
        }
    }
    MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,E->parallel.world);
    component_dot[0]=global[0];
    component_dot[1]=global[1];
}


double ala_block_vector_dot(struct All_variables *E,
                            const struct ala_block_vector *left,
                            const struct ala_block_vector *right,
                            double velocity_weight,
                            double pressure_weight)
{
    double component_dot[2];

    if(!isfinite(velocity_weight) || velocity_weight<=0.0 ||
       !isfinite(pressure_weight) || pressure_weight<=0.0)
        myerror(E,"Strict-ALA block metric weights must be positive");
    ala_block_vector_component_dot(E,left,right,component_dot);
    return velocity_weight*component_dot[0]
          +pressure_weight*component_dot[1];
}


double ala_block_vector_norm(struct All_variables *E,
                             const struct ala_block_vector *vector,
                             double velocity_weight,
                             double pressure_weight)
{
    double norm_squared;

    norm_squared=ala_block_vector_dot(E,vector,vector,
                                      velocity_weight,pressure_weight);
    if(norm_squared<0.0 && norm_squared>-1.0e-12)
        norm_squared=0.0;
    if(!isfinite(norm_squared) || norm_squared<0.0)
        myerror(E,"Invalid strict-ALA block-vector norm");
    return sqrt(norm_squared);
}


void ala_block_vector_component_norms(struct All_variables *E,
                                      const struct ala_block_vector *vector,
                                      double *velocity_norm,
                                      double *pressure_dual_norm)
{
    double component_dot[2];

    if(velocity_norm==NULL || pressure_dual_norm==NULL)
        myerror(E,"Missing strict-ALA block component norm output");
    ala_block_vector_component_dot(E,vector,vector,component_dot);
    if(!isfinite(component_dot[0]) || component_dot[0]<0.0 ||
       !isfinite(component_dot[1]) || component_dot[1]<0.0)
        myerror(E,"Invalid strict-ALA block component norm");
    *velocity_norm=sqrt(component_dot[0]);
    *pressure_dual_norm=sqrt(component_dot[1]);
}
