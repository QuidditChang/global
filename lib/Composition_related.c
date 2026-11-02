/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"
#include "parallel_related.h"
#include "composition_related.h"


static void allocate_composition_memory(struct All_variables *E);
static void compute_elemental_composition_ratio_method(struct All_variables *E);
static void compute_elemental_composition_absolute_method(struct All_variables *E);
static void init_bulk_composition(struct All_variables *E);
static void check_initial_composition(struct All_variables *E);
static void fill_composition_from_neighbors(struct All_variables *E);


/* Derived capability shared by the standalone parser and Pyre setup. */
void composition_set_capabilities(struct All_variables *E)
{
    E->composition.icompositional_rheology = E->viscosity.CDEPV;
    E->composition.on = E->control.tracer &&
        (E->composition.ichemical_buoyancy || E->composition.icompositional_rheology);
}

void composition_input(struct All_variables *E)
{
    int i;
    int m = E->parallel.me;
    input_boolean("CDEPV",&(E->viscosity.CDEPV),"off",m);

    input_boolean("chemical_buoyancy",
		  &(E->composition.ichemical_buoyancy),
		  "1,0,nomax",m);

    composition_set_capabilities(E);
    if (E->composition.on) {

        /* ibuoy_type=0 (absolute method) */
        /* ibuoy_type=1 (ratio method) */

        input_int("buoy_type",&(E->composition.ibuoy_type),"1,0,nomax",m);
        /* if (E->composition.ibuoy_type!=1) {
            fprintf(stderr,"Terror-Sorry, only ratio method allowed now\n");
            fflush(stderr);
            parallel_process_termination();
        } */

	/* Lijun add */
        E->composition.ncomp = E->trace.nflavors;

        fprintf(E->fp,"n of comp is %d\n",E->composition.icompositional_rheology);
        fprintf(stderr,"n of comp is %d\n",E->composition.ncomp);
        fprintf(E->fp,"n of flavors is %d\n",E->trace.nflavors);
        fflush(E->fp);
	/* end of add */

        if (E->composition.ibuoy_type==0)
            E->composition.ncomp = E->trace.nflavors;
        else if (E->composition.ibuoy_type==1)
            E->composition.ncomp = E->trace.nflavors - 1;
	
        E->composition.buoyancy_ratio = (double*) malloc(E->composition.ncomp
                                                         *sizeof(double));

        /* default values .... */
        for (i=0; i<E->composition.ncomp; i++)
            E->composition.buoyancy_ratio[i] = 1.0;

        input_double_vector("buoyancy_ratio", E->composition.ncomp,
                            E->composition.buoyancy_ratio,m);

    }


    /* compositional rheology */
    /* what was this about? there is a CDEPV for compositional rheology 
       
    i moved this to the top so that cdepv = on would activate tracers

    TWB  */

    /* icompositional_rheology=0 (off) */
    /* icompositional_rheology=1 (on) */
    //E->composition.icompositional_rheology = 0;


    /*
    input_int("compositional_rheology",
              &(E->composition.icompositional_rheology),"1,0,nomax",m);

    if (E->composition.icompositional_rheology==1) {
        input_double("compositional_prefactor",
                     &(E->composition.compositional_rheology_prefactor),
                     "1.0",m);
    }
    */

    return;
}



void composition_setup(struct All_variables *E)
{
    allocate_composition_memory(E);

    return;
}


void write_composition_instructions(struct All_variables *E)
{
    int k;
    double kC_min, kC_max;

    if(E->control.kC_ratio != 1.0) {
        if(!E->control.tracer || !E->composition.on ||
           E->composition.ibuoy_type != 1) {
            fprintf(stderr,
                    "Active kC requires tracers, composition, and "
                    "the ratio composition method\n");
            parallel_process_termination();
        }
        if(E->control.kC_primordial_flavor < 1 ||
           E->control.kC_primordial_flavor >= E->trace.nflavors ||
           E->control.kC_primordial_flavor - 1 >= E->composition.ncomp) {
            fprintf(stderr,
                    "kC_primordial_flavor=%d is invalid for %d tracer "
                    "flavors and %d ratio components\n",
                    E->control.kC_primordial_flavor, E->trace.nflavors,
                    E->composition.ncomp);
            parallel_process_termination();
        }
    }
    kC_min = min(1.0, E->control.kC_ratio);
    kC_max = max(1.0, E->control.kC_ratio);
    if(E->parallel.me == 0)
        fprintf(stderr,
                "Conductivity composition model: "
                "kC=1+(RC-1)*Cprim, RC=%g, primordial flavor=%d, "
                "ratio component=%d, expected kC range=[%g,%g]\n",
                E->control.kC_ratio,
                E->control.kC_primordial_flavor,
                E->control.kC_primordial_flavor - 1,
                kC_min, kC_max);

    if (E->composition.on) {

        if (E->trace.nflavors < 1) {
            fprintf(E->trace.fpt, "Tracer flavors must be greater than 1 to track composition\n");
            parallel_process_termination();
        }

        if (!E->composition.ichemical_buoyancy)
	  fprintf(E->trace.fpt,"Passive Tracers\n");
	else
	  fprintf(E->trace.fpt,"Active Tracers\n");


        if (E->composition.ibuoy_type==1)
	  fprintf(E->trace.fpt,"Ratio Method\n");
        if (E->composition.ibuoy_type==0)
	  fprintf(E->trace.fpt,"Absolute Method\n");

        for(k=0; k<E->composition.ncomp; k++) {
            fprintf(E->trace.fpt,"Buoyancy Ratio: %f\n", E->composition.buoyancy_ratio[k]);
        }

        /*
        if (E->composition.icompositional_rheology==0) {
            fprintf(E->trace.fpt,"Compositional Rheology - OFF\n");
        }
        else if (E->composition.icompositional_rheology>0) {
            fprintf(E->trace.fpt,"Compositional Rheology - ON\n");
            fprintf(E->trace.fpt,"Compositional Prefactor: %f\n",
            E->composition.compositional_rheology_prefactor);
        }
        */

        fflush(E->trace.fpt);
    }

    return;
}


/************ FILL COMPOSITION ************************/
void fill_composition(struct All_variables *E)
{

    /* ratio method */

//    if (E->composition.ibuoy_type==1) {
    if (E->composition.ibuoy_type==1 || E->composition.icompositional_rheology==1) {
        compute_elemental_composition_ratio_method(E);
    }

    /* absolute method */
    if (E->composition.ibuoy_type==0) {
        compute_elemental_composition_absolute_method(E);
    }

    /* if (E->composition.ibuoy_type!=1) {
        fprintf(E->trace.fpt,"Error(compute...)-only ratio method now\n");
        fflush(E->trace.fpt);
        exit(10);
    } */

    /* Map elemental composition to nodal points */

    map_composition_to_nodes(E);

    return;
}



static void allocate_composition_memory(struct All_variables *E)
{
    int i, j;

    for (i=0; i<E->composition.ncomp; i++) {
        E->composition.bulk_composition = (double*) malloc(E->composition.ncomp*sizeof(double));
        E->composition.initial_bulk_composition = (double*) malloc(E->composition.ncomp*sizeof(double));
        E->composition.error_fraction = (double*) malloc(E->composition.ncomp*sizeof(double));
    }

    /* for horizontal average */
    E->Have.C = (float **)malloc((E->composition.ncomp+1)*sizeof(float*));
    for (i=0; i<E->composition.ncomp; i++) {
        E->Have.C[i] = (float *)malloc((E->lmesh.noz+2)*sizeof(float));
    }


    /* allocat memory for composition fields at the nodes and elements */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        if ((E->composition.comp_el[j]=(double **)malloc((E->composition.ncomp)*sizeof(double*)))==NULL) {
            fprintf(E->trace.fpt,"AKM(allocate_composition_memory)-no memory 8987y\n");
            fflush(E->trace.fpt);
            exit(10);
        }
        if ((E->composition.comp_node[j]=(double **)malloc((E->composition.ncomp)*sizeof(double*)))==NULL) {
            fprintf(E->trace.fpt,"AKM(allocate_composition_memory)-no memory 8988y\n");
            fflush(E->trace.fpt);
            exit(10);
        }

        for (i=0; i<E->composition.ncomp; i++) {
            if ((E->composition.comp_el[j][i]=(double *)malloc((E->lmesh.nel+1)*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"AKM(allocate_composition_memory)-no memory 8989y\n");
                fflush(E->trace.fpt);
                exit(10);
            }

            if ((E->composition.comp_node[j][i]=(double *)malloc((E->lmesh.nno+1)*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"AKM(allocate_composition_memory)-no memory 983rk\n");
                fflush(E->trace.fpt);
                exit(10);
            }
        }

    }

    return;
}


void init_composition(struct All_variables *E)
{
    /* ratio method */
    if (E->composition.ibuoy_type==1 || E->composition.icompositional_rheology==1) {
        compute_elemental_composition_ratio_method(E);
    }

    /* absolute method */
    if (E->composition.ibuoy_type==0) {
        compute_elemental_composition_absolute_method(E);
    }

    /* if (E->composition.ibuoy_type!=1) {
        fprintf(E->trace.fpt,"Error(compute...)-only ratio method now\n");
        fflush(E->trace.fpt);
        exit(10);
    } */

    /* for empty elements */
    check_initial_composition(E);

    /* Map elemental composition to nodal points */
    map_composition_to_nodes(E);

    init_bulk_composition(E);
    composition_heating_initialize(E);

    return;
}


static void check_initial_composition(struct All_variables *E)
{
    /* check empty element if using ratio method */
    if (E->composition.ibuoy_type == 1) {
        if (E->trace.istat_iempty) {
            /* using the composition of neighboring elements to determine
               the initial composition of empty elements. */
            fill_composition_from_neighbors(E);
        }
    }

    return;
}



/*********** COMPUTE ELEMENTAL COMPOSITION RATIO METHOD ***/
/*                                                        */
/* This function computes the composition per element.    */
/* The concentration of material i in an element is       */
/* defined as:                                            */
/*   (# of tracers of flavor i) / (# of all tracers)      */

static void compute_elemental_composition_ratio_method(struct All_variables *E)
{
    int i, j, e, flavor, numtracers;
    int iempty = 0;


    for (j=1; j<=E->sphere.caps_per_proc; j++) {
        for (e=1; e<=E->lmesh.nel; e++) {
            numtracers = 0;
            for (flavor=0; flavor<E->trace.nflavors; flavor++)
                numtracers += E->trace.ntracer_flavor[j][flavor][e];

            /* Check for empty entries and compute ratio.  */
            /* If no tracers are in an element, skip this element, */
            /* use previous composition. */
            if (numtracers == 0) {
                iempty++;
                /* fprintf(E->trace.fpt, "No tracer in element %d!\n", e); */
                continue;
            }

            for(i=0;i<E->composition.ncomp;i++) {
                flavor = i + 1;
                E->composition.comp_el[j][i][e] =
                    E->trace.ntracer_flavor[j][flavor][e] / (double)numtracers;
            }
        }


        if (iempty) {

            if ((1.0*iempty/E->lmesh.nel)>0.80) {
                fprintf(E->trace.fpt,"WARNING(compute_elemental...)-number of tracers is REALLY LOW\n");
                fflush(E->trace.fpt);
                if (E->trace.itracer_warnings) exit(10);
            }
        }

    } /* end j */

    E->trace.istat_iempty += iempty;

    return;
}



/*********** COMPUTE ELEMENTAL ABSOLUTE METHOD *************************/
/*                                                                     */
/* This function computes the composition per element.                 */
/* The concentration of material i in an element is                    */
/* defined as:                                                         */
/*   tracer i mass * (# of tracers of flavor i) / ( volume of element) */
/*                                                                     */ 
/* where tracer i mass = tot tracer volume / tot num of tracers        */
/*                                                                     */
/* Added by DJB 01/14/10.  Known caveats: (will address at some point) */
/*     1) XXX: needs better error handling                             */
/*     2) XXX: only tested with one flavor of tracer (flavor 0)        */
/*     3) XXX: only tested by reading in a file of tracer locations    */
/*     4) XXX: the initial volume of the tracer domain is HARD-CODED!  */

static void compute_elemental_composition_absolute_method(struct All_variables *E)
{
    int i, j, e, flavor, numtracers;
    double domain_volume;
    float comp;
    float one = 1.0;

    /* This needs to be `manually' changed for the total volume */
    /*  occupied by your tracers */
    domain_volume = 1e-2;

    for (j=1; j<=E->sphere.caps_per_proc; j++) {
        for (e=1; e<=E->lmesh.nel; e++) {
            numtracers = 0;
            for (flavor=0; flavor<E->trace.nflavors; flavor++)
                numtracers += E->trace.ntracer_flavor[j][flavor][e];

            /* Check for empty entries */
            /* If no tracers are in an element, comp = 0.0 (i.e. is ambient) */
            if (numtracers == 0) {
                for(i=0;i<E->composition.ncomp;i++) {
                    E->composition.comp_el[j][i][e] = 0.0;
                }
                continue;
            }

            /* Composition is proportional to (local) tracer density */
            for(i=0;i<E->composition.ncomp;i++) {
                flavor = i;
                comp =
                    E->trace.ntracer_flavor[j][flavor][e] / E->eco[j][e].area
                    * domain_volume / E->trace.number_of_tracers;

                /* truncate composition at 1.0 */
                /* This violates mass conservation but prevents unphysical C */
                /* XXX: make truncation a switch for the user to specify */
                E->composition.comp_el[j][i][e] = min(comp,one);

            }
        }
    }

    return;
}



/********** MAP COMPOSITION TO NODES ****************/
/*                                                  */


void map_composition_to_nodes(struct All_variables *E)
{
    double *tmp[NCS];
    int i, n, kk;
    int nelem, nodenum;
    int j;


    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        /* first, initialize node array */
        for(i=0;i<E->composition.ncomp;i++) {
            for (kk=1;kk<=E->lmesh.nno;kk++)
                E->composition.comp_node[j][i][kk]=0.0;
        }

        /* Loop through all elements */
        for (nelem=1;nelem<=E->lmesh.nel;nelem++) {

            /* for each element, loop through element nodes */

            /* weight composition */

            for (nodenum=1;nodenum<=8;nodenum++) {
                n = E->ien[j][nelem].node[nodenum];
                for(i=0;i<E->composition.ncomp;i++) {

                    E->composition.comp_node[j][i][n] +=
                        E->composition.comp_el[j][i][nelem]*
                        E->TWW[E->mesh.levmax][j][nelem].node[nodenum];
                }
            }

        } /* end nelem */
    } /* end j */

    for(i=0;i<E->composition.ncomp;i++) {
        for (j=1;j<=E->sphere.caps_per_proc;j++)
            tmp[j] = E->composition.comp_node[j][i];

        (E->exchange_node_d)(E,tmp,E->mesh.levmax);
    }

    /* Divide by nodal volume */
    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for(i=0;i<E->composition.ncomp;i++)
            for (kk=1;kk<=E->lmesh.nno;kk++)
                E->composition.comp_node[j][i][kk] *= E->MASS[E->mesh.levmax][j][kk];

        /* testing */
        /**
        for(i=0;i<E->composition.ncomp;i++)
            for (kk=1;kk<=E->lmesh.nel;kk++) {
                fprintf(E->trace.fpt,"%d %f\n",kk,E->composition.comp_el[j][i][kk]);
            }

        for(i=0;i<E->composition.ncomp;i++)
            for (kk=1;kk<=E->lmesh.nno;kk++) {
                fprintf(E->trace.fpt,"%d %f %f\n",kk,E->sx[j][3][kk],E->composition.comp_node[j][i][kk]);
            }
        fflush(E->trace.fpt);
        /**/

    } /* end j */

    return;
}


/****************************************************************/

static void fill_composition_from_neighbors(struct All_variables *E)
{
    int i, j, k, e, ee, n, flavor, numtracers, count;
    double *sum;
    const int n_nghbrs = 4;
    int nghbrs[n_nghbrs];
    int *is_empty;

    fprintf(E->trace.fpt,"WARNING(check_initial_composition)-number of tracers is low, %d elements contain no tracer initially\n", E->trace.istat_iempty);

    fprintf(E->trace.fpt,"Using neighboring elements for initial composition...\n");

    /* index shift for neighboring elements in horizontal direction */
    nghbrs[0] = E->lmesh.elz;
    nghbrs[1] = -E->lmesh.elz;
    nghbrs[2] = E->lmesh.elz * E->lmesh.elx;
    nghbrs[3] = -E->lmesh.elz * E->lmesh.elx;

    is_empty = (int *)calloc(E->lmesh.nel+1, sizeof(int));
    sum = (double *)malloc(E->composition.ncomp * sizeof(double));

    for (j=1; j<=E->sphere.caps_per_proc; j++) {
        /* which element is empty? */
        for (e=1; e<=E->lmesh.nel; e++) {
            numtracers = 0;
            for (flavor=0; flavor<E->trace.nflavors; flavor++)
                numtracers += E->trace.ntracer_flavor[j][flavor][e];

            if (numtracers == 0)
                is_empty[e] = 1;
        }

        /* using the average comp_el from neighboring elements */
        for (e=1; e<=E->lmesh.nel; e++) {
            if(is_empty[e]) {
                count = 0;
                for (i=0; i<E->composition.ncomp; i++)
                    sum[i] = 0.0;

                for(n=0; n<n_nghbrs; n++) {
                    ee = e + nghbrs[n];
                    /* is ee a valid element number and the elemnt is not empty? */
                    if((ee>0) && (ee<=E->lmesh.nel) && (!is_empty[ee])) {
                        count++;
                        for (i=0; i<E->composition.ncomp; i++)
                            sum[i] += E->composition.comp_el[j][i][ee];
                    }
                }

                if(count == 0) {
                    fprintf(E->trace.fpt,"Error(fill_composition_from_neighbors)-all neighboring elements are empty\n");
		    fprintf(E->trace.fpt,"rad:%f,%f,%f,%f,%f,%f",E->sx[j][3][E->ien[j][e].node[1]],E->sx[j][3][E->ien[j][e].node[8]],E->sx[j][1][E->ien[j][e].node[1]],E->sx[j][1][E->ien[j][e].node[3]],E->sx[j][2][E->ien[j][e].node[1]],E->sx[j][2][E->ien[j][e].node[3]]);
                    fflush(E->trace.fpt);
                    exit(10);
                }

                for (i=0; i<E->composition.ncomp; i++)
                    E->composition.comp_el[j][i][e] = sum[i] / count;
            }
        }
    }

    free(is_empty);
    free(sum);

    fprintf(E->trace.fpt,"Done.\n");
    fflush(E->trace.fpt);
    return;
}


/*********** GET BULK COMPOSITION *******************************/

static void init_bulk_composition(struct All_variables *E)
{

    double return_bulk_value_d();
    double volume;
    double *tmp[NCS];
    int i, m;
    const int ival=0;


    for (i=0; i<E->composition.ncomp; i++) {

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            tmp[m] = E->composition.comp_node[m][i];

        /* ival=0 returns integral not average */
        volume = return_bulk_value_d(E,tmp,ival);

        E->composition.bulk_composition[i] = volume;
        E->composition.initial_bulk_composition[i] = volume;
    }

    return;
}


void get_bulk_composition(struct All_variables *E)
{

    double return_bulk_value_d();
    double volume;
    double *tmp[NCS];
    int i, m;
    const int ival = 0;

    for (i=0; i<E->composition.ncomp; i++) {

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            tmp[m] = E->composition.comp_node[m][i];

        /* ival=0 returns integral not average */
        volume = return_bulk_value_d(E,tmp,ival);

        E->composition.bulk_composition[i] = volume;

        E->composition.error_fraction[i] = (volume - E->composition.initial_bulk_composition[i]) / E->composition.initial_bulk_composition[i];
    }

    return;
}

/* RH is a material contrast, normalized once with the actual source quadrature.
 * State is checkpointed; no normalization is performed on advancing composition. */
static void rh_error(struct All_variables *E,const char *reason)
{
    fprintf(stderr,"RH_ERROR rank=%d: %s\n",E->parallel.me,reason);
    MPI_Abort(E->parallel.world,72);
}

static void rh_integrals(struct All_variables *E,double total[4])
{
    int m,e,z,g,f,count; double local[4]={0,0,0,0};
    struct Shape_function GN; struct Shape_function_dx dx;
    struct Shape_function_dA omega; double rtf[4][9];
    void get_global_shape_fn();
    for(m=1;m<=E->sphere.caps_per_proc;m++)for(e=1;e<=E->lmesh.nel;e++) {
        double v=0,c=E->composition.comp_el[m][23][e],rho;
        if(!isfinite(c)||c<0||c>1)rh_error(E,"invalid primordial fraction");
        get_global_shape_fn(E,e,&GN,&dx,&omega,0,1,rtf,E->mesh.levmax,m);
        for(g=1;g<=8;g++)v+=omega.vpt[g]*g_point[g].weight[2];
        z=(e-1)%E->lmesh.elz+1;
        rho=.5*(E->refstate.rho[z]+E->refstate.rho[z+1]);
        if(!isfinite(v*rho)||v<=0||rho<=0)rh_error(E,"invalid source integration weight");
        local[0]+=rho*v; local[1]+=rho*v*c;
        count=0;for(f=0;f<E->trace.nflavors;f++)count+=E->trace.ntracer_flavor[m][f][e];
        if(!count)local[2]+=rho*v*c;
        local[3]+=v*c;
    }
    MPI_Allreduce(local,total,4,MPI_DOUBLE,MPI_SUM,E->parallel.world);
}

void composition_heating_initialize(struct All_variables *E)
{
    double t[4],*s=E->control.rh_state;
    if(!isfinite(E->control.RH)||E->control.RH<1)rh_error(E,"RH must be finite and >= 1");
    if(E->control.RH==1)return;
    if(!E->pices.enabled || !E->pices.checkpoint || !E->composition.on ||
       E->composition.ibuoy_type!=1 || E->composition.ncomp<24 ||
       E->control.kC_primordial_flavor!=24 || E->control.tracer_enriched)
        rh_error(E,"RH requires PICES checkpointing and ratio flavor 24; legacy enrichment is incompatible");
    if(!isfinite(E->control.Q0)||E->control.Q0<0)rh_error(E,"RH requires finite nonnegative Q0");
    if(E->control.restart)return; /* restored from the authenticated sidecar */
    if(E->control.rh_initialized)return;
    rh_integrals(E,t);
    s[0]=t[1]/t[0];s[1]=E->control.Q0/(1+(E->control.RH-1)*s[0]);
    s[2]=E->control.RH*s[1];s[3]=E->control.Q0*t[0];
    E->control.rh_initialized=1;
    composition_heating_report(E);
}

double composition_heating(struct All_variables *E,int m,int e)
{
    double c;
    if(E->control.RH==1)return E->control.Q0;
    if(!E->control.rh_initialized)rh_error(E,"missing fixed heating coefficients");
    c=E->composition.comp_el[m][23][e];
    if(!isfinite(c)||c<0||c>1)rh_error(E,"invalid primordial fraction");
    return (1-c)*E->control.rh_state[1]+c*E->control.rh_state[2];
}

void composition_heating_report(struct All_variables *E)
{
    double t[4],*s=E->control.rh_state,scale,b,p;
    if(E->control.RH==1)return;
    rh_integrals(E,t);
    scale=E->data.k0*E->data.ref_temperature*E->data.radius_km*1000.;
    b=(t[0]-t[1])*s[1];p=t[1]*s[2];
    if(E->parallel.me==0)fprintf(E->fp,
      "RH_HEATING step=%d RH=%.17g f0=%.17g Qbg=%.17g Qprim=%.17g f_current=%.17g primordial_volume_nd=%.17g background_W=%.17g primordial_W=%.17g total_W=%.17g baseline_W=%.17g drift_W=%.17g empty_primordial_W=%.17g scope=whole_shell state=current_composition\n",
      E->monitor.solution_cycles,(double)E->control.RH,s[0],s[1],s[2],t[1]/t[0],t[3],b*scale,p*scale,(b+p)*scale,s[3]*scale,(b+p-s[3])*scale,t[2]*s[2]*scale);
}
