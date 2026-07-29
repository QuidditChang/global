/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include "global_defs.h"

#include "parallel_related.h"
#include "parsing.h"
#include "phase_change.h"

static void apply_one_phase(struct All_variables *E, double **buoy,
                            int phase_index);
static void calc_phase_change(struct All_variables *E,
                              int phase_index);
static void debug_phase_change(struct All_variables *E, int phase_index);
static void validate_phase_parameters(struct All_variables *E);


void phase_change_allocate(struct All_variables *E)
{
  int j, phase_index;
  int nno  = E->lmesh.nno;
  int nsf  = E->lmesh.nsf;

  validate_phase_parameters(E);

  for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
    for(j=1; j<=E->sphere.caps_per_proc; j++) {
      E->phase_B[phase_index][j] =
          (float *) malloc((nno+1)*sizeof(float));
      E->phase_boundary[phase_index][j] =
          (float *) malloc((nsf+1)*sizeof(float));
    }

  return;
}


void phase_change_input(struct All_variables *E)
{
  int m = E->parallel.me;
  int phase_index;
  float depth[PHASE_TRANSITIONS];
  float density_jump[PHASE_TRANSITIONS];
  float Ra[PHASE_TRANSITIONS];
  float width[PHASE_TRANSITIONS];
  float clapeyron[PHASE_TRANSITIONS];
  float transT[PHASE_TRANSITIONS];

  input_float_vector("phase_depth", PHASE_TRANSITIONS, depth, m);
  input_float_vector("phase_delta_rho", PHASE_TRANSITIONS, density_jump, m);
  input_float_vector("phase_Ra", PHASE_TRANSITIONS, Ra, m);
  input_float_vector("phase_width", PHASE_TRANSITIONS, width, m);
  input_float_vector("phase_clapeyron", PHASE_TRANSITIONS, clapeyron, m);
  input_float_vector("phase_transT", PHASE_TRANSITIONS, transT, m);

  for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++) {
    E->control.phase[phase_index].depth = depth[phase_index];
    E->control.phase[phase_index].density_jump = density_jump[phase_index];
    E->control.phase[phase_index].Ra = Ra[phase_index];
    E->control.phase[phase_index].clapeyron = clapeyron[phase_index];
    E->control.phase[phase_index].transT = transT[phase_index];
    E->control.phase[phase_index].inv_width =
        (width[phase_index] == 0.0)? 0.0 : 1.0/width[phase_index];
  }

  return;
}


void phase_change_apply(struct All_variables *E, double **buoy)
{
  int phase_index;

  for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++)
    if(E->control.phase[phase_index].Ra != 0.0)
      apply_one_phase(E, buoy, phase_index);

  return;
}


float phase_change_reference_fraction(struct All_variables *E,
                                      int phase_index, int cap, int node)
{
  int nz;
  float e_pressure, dz;
  struct Phase_transition *phase = &E->control.phase[phase_index];

  nz = ((node-1) % E->lmesh.noz) + 1;
  dz = (E->sphere.ro-E->sx[cap][3][node]) - phase->depth;
  e_pressure = dz * E->refstate.rho[nz] * E->refstate.gravity[nz]
      - phase->clapeyron
      * (E->refstate.temperature[nz] - phase->transT);

  return 0.5 * (1.0 + tanh(phase->inv_width * e_pressure));
}


static void apply_one_phase(struct All_variables *E, double **buoy,
                            int phase_index)
{
  int m, i;
  float Xref;
  struct Phase_transition *phase = &E->control.phase[phase_index];
  float **B = E->phase_B[phase_index];

  calc_phase_change(E, phase_index);
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++) {
      Xref = phase_change_reference_fraction(E, phase_index, m, i);
      buoy[m][i] -= phase->Ra * (B[m][i] - Xref);
    }

  if (E->control.verbose) {
    fprintf(E->fp_out,
            "phase[%d] depth=%f depth_km=%f delta_rho=%f Ra=%f "
            "clapeyron=%f transT=%f inv_width=%f\n",
            phase_index, phase->depth, phase->depth*E->data.radius_km,
            phase->density_jump, phase->Ra, phase->clapeyron,
            phase->transT, phase->inv_width);
    debug_phase_change(E, phase_index);
    fflush(E->fp_out);
  }

  return;
}


static void calc_phase_change(struct All_variables *E,
                              int phase_index)
{
  int i,j,k,n,ns,m,nz;
  float e_pressure,pt5,one,dz;
  struct Phase_transition *phase = &E->control.phase[phase_index];
  float **B = E->phase_B[phase_index];
  float **B_b = E->phase_boundary[phase_index];

  pt5 = 0.5;
  one = 1.0;

  for(m=1;m<=E->sphere.caps_per_proc;m++)     {
    /* compute phase function B, the concentration of the high pressure
     * phase. B is between 0 and 1. */
    for(i=1;i<=E->lmesh.nno;i++)  {
        nz = ((i-1) % E->lmesh.noz) + 1;
        dz = (E->sphere.ro-E->sx[m][3][i]) - phase->depth;
        /*XXX: dz*rho[nz]*g[nz] is only a approximation for the reduced
         * pressure, a more accurate formula is:
         *   integral(rho(z)*g(z)*dz) from depth_ph to current depth   */
        e_pressure = dz * E->refstate.rho[nz] * E->refstate.gravity[nz]
            - phase->clapeyron * (E->T[m][i] - phase->transT);

        B[m][i] = pt5 * (one + tanh(phase->inv_width * e_pressure));
    }

    /* compute the phase boundary, defined as the depth where B==0.5 */
    ns = 0;
    for (k=1;k<=E->lmesh.noy;k++)
      for (j=1;j<=E->lmesh.nox;j++)  {
        ns++;
        B_b[m][ns]=0.0;
        for (i=1;i<E->lmesh.noz;i++)   {
          n = (k-1)*E->lmesh.noz*E->lmesh.nox + (j-1)*E->lmesh.noz + i;
          if (B[m][n]>=pt5 && B[m][n+1]<=pt5)
            B_b[m][ns]=(E->sx[m][3][n+1]-E->sx[m][3][n])*(pt5-B[m][n])/(B[m][n+1]-B[m][n])+E->sx[m][3][n];
	}
      }
  }

  return;
}


static void debug_phase_change(struct All_variables *E, int phase_index)
{
  int m, j;
  float Xref, deltaX;
  struct Phase_transition *phase = &E->control.phase[phase_index];
  float **B = E->phase_B[phase_index];

  fprintf(E->fp_out,
          "output_phase_change phase_index=%d phase_depth=%f depth_km=%f\n",
          phase_index, phase->depth, phase->depth*E->data.radius_km);
  for(m=1;m<=E->sphere.caps_per_proc;m++)        {
    fprintf(E->fp_out,"for cap %d\n",E->sphere.capid[m]);
    for (j=1;j<=E->lmesh.nno;j++) {
      Xref = phase_change_reference_fraction(E, phase_index, m, j);
      deltaX = B[m][j] - Xref;
      fprintf(E->fp_out,
              "phase_depth=%.6e node=%06d Z=%.6e T=%.6e "
              "B=%.6e Xref=%.6e deltaX=%.6e "
              "density_jump_contribution=%.6e "
              "phase_buoyancy_contribution=%.6e\n",
              phase->depth, j, E->sx[m][3][j], E->T[m][j], B[m][j],
              Xref, deltaX, phase->density_jump*deltaX,
              -phase->Ra*deltaX);
    }
  }
  fflush(E->fp_out);

  return;
}


static void validate_phase_parameters(struct All_variables *E)
{
  int phase_index;
  struct Phase_transition *phase;
  double density_scale, expected_Ra, relative_error;

  density_scale = E->data.density * E->data.therm_exp
      * E->data.ref_temperature;

  for(phase_index=0; phase_index<PHASE_TRANSITIONS; phase_index++) {
    phase = &E->control.phase[phase_index];

    if(phase->depth < 0.0 || phase->depth > E->sphere.ro-E->sphere.ri) {
      fprintf(stderr, "phase[%d] depth=%g is outside the model shell\n",
              phase_index, phase->depth);
      parallel_process_termination();
    }

    if(phase->Ra != 0.0 && phase->inv_width == 0.0) {
      fprintf(stderr, "phase[%d] has nonzero Ra but zero width\n", phase_index);
      parallel_process_termination();
    }

    if(phase->density_jump != 0.0 && phase->Ra != 0.0 &&
       density_scale != 0.0 && E->control.Atemp != 0.0) {
      expected_Ra = E->control.Atemp * phase->density_jump / density_scale;
      relative_error = fabs(phase->Ra-expected_Ra) / fabs(expected_Ra);
      if(relative_error > 0.01 && E->parallel.me == 0)
        fprintf(stderr,
                "phase[%d] Ra=%g differs from delta_rho-derived Ra=%g "
                "by %.2f percent\n",
                phase_index, phase->Ra, expected_Ra, 100.0*relative_error);
    }
  }

  return;
}
