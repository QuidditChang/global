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
#include "element_definitions.h"
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
          (float *) calloc(nno+1, sizeof(float));
      E->phase_boundary[phase_index][j] =
          (float *) calloc(nsf+1, sizeof(float));
    }

  return;
}


void phase_change_input(struct All_variables *E)
{
  int m = E->parallel.me;
  int phase_index;
  float depth[LEGACY_PHASE_TRANSITIONS];
  float density_jump[LEGACY_PHASE_TRANSITIONS];
  float entropy_jump[LEGACY_PHASE_TRANSITIONS];
  float width[LEGACY_PHASE_TRANSITIONS];
  float clapeyron[LEGACY_PHASE_TRANSITIONS];
  float transT[LEGACY_PHASE_TRANSITIONS];
  float pPv_depth, pPv_transT, pPv_clapeyron, pPv_width;
  float pPv_density, pPv_entropy;

  input_float_vector("phase_depth", LEGACY_PHASE_TRANSITIONS, depth, m);
  input_float_vector("phase_delta_rho", LEGACY_PHASE_TRANSITIONS, density_jump, m);
  if(!input_float_vector("phase_delta_s", LEGACY_PHASE_TRANSITIONS, entropy_jump, m)) {
    fprintf(stderr, "phase_delta_s is required for the single-entropy phase model\n");
    parallel_process_termination();
  }
  input_float_vector("phase_width", LEGACY_PHASE_TRANSITIONS, width, m);
  input_float_vector("phase_clapeyron", LEGACY_PHASE_TRANSITIONS, clapeyron, m);
  input_float_vector("phase_transT", LEGACY_PHASE_TRANSITIONS, transT, m);

  for(phase_index=0; phase_index<LEGACY_PHASE_TRANSITIONS; phase_index++) {
    E->control.phase[phase_index].depth = depth[phase_index];
    E->control.phase[phase_index].density_jump = density_jump[phase_index];
    E->control.phase[phase_index].entropy_jump = entropy_jump[phase_index];
    E->control.phase[phase_index].Ra = 0.0;
    E->control.phase[phase_index].clapeyron = clapeyron[phase_index];
    E->control.phase[phase_index].transT = transT[phase_index];
    E->control.phase[phase_index].inv_width =
        (width[phase_index] == 0.0)? 0.0 : 1.0/width[phase_index];
  }

  input_boolean("pPv_enabled", &E->control.pPv_enabled, "off", m);
  input_boolean("pPv_latent_enabled", &E->control.pPv_latent_enabled, "off", m);
  input_float("pPv_k_ratio", &E->control.pPv_k_ratio, "1.0", m);
  input_float("pPv_eta_ratio", &E->control.pPv_eta_ratio, "1.0", m);
  input_float("pPv_depth", &pPv_depth, "0.0", m);
  input_float("pPv_transT", &pPv_transT, "0.0", m);
  input_float("pPv_clapeyron", &pPv_clapeyron, "0.0", m);
  input_float("pPv_width", &pPv_width, "0.0", m);
  input_float("pPv_delta_rho", &pPv_density, "0.0", m);
  input_float("pPv_delta_s", &pPv_entropy, "0.0", m);
  phase_change_set_pPv(E, pPv_depth, pPv_transT, pPv_clapeyron,
                       pPv_width, pPv_density, pPv_entropy);

  return;
}


void phase_change_set_pPv(struct All_variables *E, float depth, float transT,
                          float clapeyron, float width, float density_jump,
                          float entropy_jump)
{
  struct Phase_transition *phase = &E->control.phase[pPv_PHASE_INDEX];
  phase->depth = depth;
  phase->transT = transT;
  phase->clapeyron = clapeyron;
  phase->inv_width = width == 0.0 ? 0.0 : 1.0/width;
  phase->density_jump = density_jump;
  phase->entropy_jump = entropy_jump;
  phase->Ra = 0.0;
}


int phase_change_count(const struct All_variables *E)
{
  return E->control.pPv_enabled ? PHASE_TRANSITIONS : LEGACY_PHASE_TRANSITIONS;
}


double phase_change_entropy(const struct All_variables *E, int phase_index)
{
  if(phase_index == pPv_PHASE_INDEX &&
     (!E->control.pPv_enabled || !E->control.pPv_latent_enabled))
    return 0.0;
  return E->control.phase[phase_index].entropy_jump;
}


int phase_change_thermal_active(const struct All_variables *E, int phase_index)
{
  if(phase_index < LEGACY_PHASE_TRANSITIONS)
    return 1; /* Preserve the existing three-phase substep constraints. */
  return E->control.pPv_enabled &&
      (E->control.pPv_latent_enabled || E->control.pPv_k_ratio != 1.0);
}


void phase_change_apply(struct All_variables *E, double **buoy)
{
  int phase_index;

  for(phase_index=0; phase_index<phase_change_count(E); phase_index++)
    if(E->control.phase[phase_index].Ra != 0.0)
      apply_one_phase(E, buoy, phase_index);
    else
      calc_phase_change(E, phase_index);

  return;
}


void phase_change_state(const struct Phase_transition *phase,
                        double depth, double temperature,
                        double rho, double gravity,
                        double d_rho_g_dr,
                        double *q, double *fraction,
                        double *dX_dT, double *dX_dr_pressure)
{
  double argument, localization, rho_g;

  rho_g = rho * gravity;
  *q = (depth - phase->depth) * rho_g
      - phase->clapeyron * (temperature - phase->transT);
  argument = phase->inv_width * *q;
  *fraction = 0.5 * (1.0 + tanh(argument));
  localization = 2.0 * phase->inv_width * *fraction * (1.0 - *fraction);
  *dX_dT = -localization * phase->clapeyron;
  *dX_dr_pressure = localization
      * (-rho_g + (depth - phase->depth) * d_rho_g_dr);
}


float phase_change_reference_fraction(const struct All_variables *E,
                                      int phase_index, int cap, int node)
{
  int nz;
  double q, fraction, dX_dT, dX_dr_pressure;
  const struct Phase_transition *phase = &E->control.phase[phase_index];

  if(phase_index == pPv_PHASE_INDEX && !E->control.pPv_enabled)
    return 0.0;

  nz = ((node-1) % E->lmesh.noz) + 1;
  phase_change_state(phase, E->sphere.ro-E->sx[cap][3][node],
                     E->refstate.Tref[nz], E->refstate.rho[nz],
                     E->refstate.gravity[nz], 0.0, &q, &fraction,
                     &dX_dT, &dX_dr_pressure);
  return fraction;
}


float phase_change_fraction_at_temperature(const struct All_variables *E,
                                           int phase_index, int cap, int node,
                                           double temperature)
{
  int nz;
  double q, fraction, dX_dT, dX_dr_pressure;
  const struct Phase_transition *phase = &E->control.phase[phase_index];

  if(phase_index == pPv_PHASE_INDEX && !E->control.pPv_enabled)
    return 0.0;

  nz = ((node-1) % E->lmesh.noz) + 1;
  phase_change_state(phase, E->sphere.ro-E->sx[cap][3][node], temperature,
                     E->refstate.rho[nz], E->refstate.gravity[nz], 0.0,
                     &q, &fraction, &dX_dT, &dX_dr_pressure);
  return fraction;
}


double phase_change_density_fraction(const struct All_variables *E,
                                      int phase_index, int cap, int node,
                                      double fraction)
{
  if(phase_index == pPv_PHASE_INDEX)
    return E->control.pPv_enabled ? fraction : 0.0;
  return fraction - phase_change_reference_fraction(E, phase_index, cap, node);
}


double pPv_fraction_at_state(const struct All_variables *E, double depth,
                             double temperature, double rho, double gravity)
{
  double q, fraction, dT, dr;
  if(!E->control.pPv_enabled)
    return 0.0;
  phase_change_state(&E->control.phase[pPv_PHASE_INDEX], depth, temperature,
                     rho, gravity, 0.0, &q, &fraction, &dT, &dr);
  return fraction;
}


double pPv_fraction_at_node(const struct All_variables *E, int cap, int node,
                            double temperature)
{
  int nz;
  if(!E->control.pPv_enabled)
    return 0.0;
  nz = ((node-1) % E->lmesh.noz) + 1;
  return pPv_fraction_at_state(E, E->sphere.ro-E->sx[cap][3][node],
      temperature, E->refstate.rho[nz], E->refstate.gravity[nz]);
}


double pPv_fraction_at_gp(const struct All_variables *E, int cap, int element,
                          int gp, double temperature)
{
  int a, j, node, nz;
  double xyz[3] = {0.0, 0.0, 0.0};
  double rho_g = 0.0, weight, radius;
  if(!E->control.pPv_enabled)
    return 0.0;
  for(a=1; a<=8; a++) {
    node = E->ien[cap][element].node[a];
    nz = ((node-1) % E->lmesh.noz) + 1;
    weight = E->N.vpt[GNVINDEX(a,gp)];
    for(j=0; j<3; j++)
      xyz[j] += weight * E->x[cap][j+1][node];
    rho_g += weight * E->refstate.rho[nz] * E->refstate.gravity[nz];
  }
  radius = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]+xyz[2]*xyz[2]);
  return pPv_fraction_at_state(E, E->sphere.ro-radius, temperature, 1.0, rho_g);
}


double pPv_conductivity_factor(const struct All_variables *E, double fraction)
{
  if(!E->control.pPv_enabled || E->control.pPv_k_ratio == 1.0)
    return 1.0;
  return 1.0 + (E->control.pPv_k_ratio-1.0)*fraction;
}


double pPv_viscosity_factor(const struct All_variables *E, double fraction)
{
  if(!E->control.pPv_enabled || E->control.pPv_eta_ratio == 1.0)
    return 1.0;
  return pow(E->control.pPv_eta_ratio, fraction);
}


static void apply_one_phase(struct All_variables *E, double **buoy,
                            int phase_index)
{
  int m, i;
  struct Phase_transition *phase = &E->control.phase[phase_index];
  float **B = E->phase_B[phase_index];

  calc_phase_change(E, phase_index);
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++) {
      if(phase_index < LEGACY_PHASE_TRANSITIONS) {
        /* Preserve the original float arithmetic for the legacy solver path. */
        float Xref = phase_change_reference_fraction(E, phase_index, m, i);
        buoy[m][i] -= phase->Ra * (B[m][i] - Xref);
      }
      else
        buoy[m][i] -= phase->Ra * phase_change_density_fraction(
            E, phase_index, m, i, B[m][i]);
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
  float pt5;
  double q, fraction, dX_dT, dX_dr_pressure;
  struct Phase_transition *phase = &E->control.phase[phase_index];
  float **B = E->phase_B[phase_index];
  float **B_b = E->phase_boundary[phase_index];

  pt5 = 0.5;

  for(m=1;m<=E->sphere.caps_per_proc;m++)     {
    /* compute phase function B, the concentration of the high pressure
     * phase. B is between 0 and 1. */
    for(i=1;i<=E->lmesh.nno;i++)  {
        nz = ((i-1) % E->lmesh.noz) + 1;
        phase_change_state(phase, E->sphere.ro-E->sx[m][3][i], E->T[m][i],
                           E->refstate.rho[nz], E->refstate.gravity[nz], 0.0,
                           &q, &fraction, &dX_dT, &dX_dr_pressure);
        B[m][i] = fraction;
    }

    /* The legacy cache represents one crossing only. pPv may have zero or
     * multiple crossings; its full fraction field is the authoritative output. */
    if(phase_index == pPv_PHASE_INDEX)
      continue;

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
      deltaX = phase_change_density_fraction(E, phase_index, m, j, B[m][j]);
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
  double density_scale, expected_Ra, legacy_Ra, relative_error;
  double pressure_scale_gpa, transition_temperature_k, clapeyron_mpa_k;
  double width_gpa;

  density_scale = E->data.rho0 * E->data.alpha0
      * E->data.ref_temperature;
  pressure_scale_gpa = E->data.rho0 * E->data.g0
      * E->data.radius_km * 1.0e-6;

  if(!isfinite(density_scale) || density_scale <= 0.0 ||
     !isfinite(E->control.Atemp)) {
    fprintf(stderr, "cannot derive phase_Ra from invalid reference scales\n");
    parallel_process_termination();
  }

  if(E->control.pPv_enabled &&
     (!isfinite(E->control.pPv_k_ratio) || E->control.pPv_k_ratio <= 0.0 ||
      !isfinite(E->control.pPv_eta_ratio) || E->control.pPv_eta_ratio <= 0.0)) {
    fprintf(stderr, "pPv conductivity and viscosity ratios must be positive and finite\n");
    parallel_process_termination();
  }

  for(phase_index=0; phase_index<phase_change_count(E); phase_index++) {
    phase = &E->control.phase[phase_index];
    legacy_Ra = phase->Ra;
    expected_Ra = E->control.Atemp * phase->density_jump / density_scale;
    phase->Ra = expected_Ra;

    if(!isfinite(phase->depth) || phase->depth < 0.0 ||
       phase->depth > E->sphere.ro-E->sphere.ri) {
      fprintf(stderr, "phase[%d] depth=%g is outside the model shell\n",
              phase_index, phase->depth);
      parallel_process_termination();
    }

    if(!isfinite(phase->density_jump) || !isfinite(phase->clapeyron) ||
       !isfinite(phase->transT) || !isfinite(phase->inv_width) ||
       phase->inv_width < 0.0 ||
       ((expected_Ra != 0.0 || phase->entropy_jump != 0.0 ||
         phase_index == pPv_PHASE_INDEX) && phase->inv_width <= 0.0)) {
      fprintf(stderr, "phase[%d] has nonfinite parameters or nonpositive active width\n", phase_index);
      parallel_process_termination();
    }

    if(!isfinite(phase->entropy_jump)) {
      fprintf(stderr, "phase[%d] has non-finite entropy jump\n", phase_index);
      parallel_process_termination();
    }
    if(phase_index < LEGACY_PHASE_TRANSITIONS &&
       phase->density_jump != 0.0 && phase->entropy_jump == 0.0) {
      fprintf(stderr, "phase[%d] has density contrast but no entropy jump\n",
              phase_index);
      parallel_process_termination();
    }

    if(phase_index == pPv_PHASE_INDEX &&
       (phase->density_jump < 0.0 ||
        E->data.Ttop + phase->transT*E->data.ref_temperature <= 0.0 ||
        (E->control.pPv_latent_enabled && phase->density_jump > 0.0 &&
         phase->clapeyron != 0.0 &&
         (phase->entropy_jump == 0.0 || phase->clapeyron*phase->entropy_jump > 0.0)))) {
      fprintf(stderr, "pPv requires nonnegative density increment, positive anchor temperature, and consistent latent entropy\n");
      parallel_process_termination();
    }

    if(legacy_Ra != 0.0) {
      relative_error = fabs(legacy_Ra-expected_Ra)
          / fmax(fabs(expected_Ra), 1.0);
      if(relative_error > 0.01 && E->parallel.me == 0)
        fprintf(stderr,
                "phase[%d] legacy phase_Ra=%g ignored; "
                "delta_rho-derived phase_Ra=%g\n",
                phase_index, legacy_Ra, expected_Ra);
    }

    if(E->parallel.me == 0) {
      transition_temperature_k = E->data.Ttop
          + phase->transT * E->data.ref_temperature;
      clapeyron_mpa_k = phase->clapeyron * pressure_scale_gpa * 1.0e3
          / E->data.ref_temperature;
      width_gpa = (phase->inv_width == 0.0)
          ? 0.0 : pressure_scale_gpa / phase->inv_width;
      fprintf(stderr,
              "phase[%d] z0_km=%g T0_K=%g Gamma_MPa_K=%g "
              "width_GPa=%g delta_rho_kg_m3=%g delta_s_J_kgK=%g "
              "derived_phase_Ra=%g\n",
              phase_index, phase->depth * E->data.radius_km,
              transition_temperature_k, clapeyron_mpa_k, width_gpa,
              phase->density_jump, phase->entropy_jump * E->data.Cp0,
              phase->Ra);
    }
  }

  if(E->parallel.me == 0 && E->control.pPv_enabled)
    fprintf(stderr, "pPv enabled: direct density increment; latent=%d k_ratio=%g eta_ratio=%g\n",
            E->control.pPv_latent_enabled, E->control.pPv_k_ratio,
            E->control.pPv_eta_ratio);

  return;
}
