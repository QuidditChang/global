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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <math.h>
#include <string.h>

#include "global_defs.h"
#include "material_properties.h"
#include "parallel_related.h"

static void read_refstate(struct All_variables *E);
static void read_ala_beta_intervals(struct All_variables *E);
static int read_refstate_data_line(FILE *fp, char *buffer, int length);
static void adams_williamson_eos(struct All_variables *E);
static void initialize_ala_beta(struct All_variables *E);

int layers_r(struct All_variables *,float);

void mat_prop_allocate(struct All_variables *E)
{
    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;

    /* reference profile of density */
    E->refstate.rho = (double *) malloc((noz+1)*sizeof(double));

    /* Strict ALA nodal input and its one authoritative element mapping. */
    E->refstate.ala_beta =
        (double *) malloc((E->lmesh.elz+1)*sizeof(double));
    E->refstate.ala_beta_supplied =
        (double *) malloc((E->lmesh.elz+1)*sizeof(double));
    E->refstate.ala_beta_density =
        (double *) malloc((E->lmesh.elz+1)*sizeof(double));
    E->refstate.ala_beta_interval =
        (double *) calloc(E->lmesh.elz+1, sizeof(double));
    E->refstate.beta_ala = (double *) calloc(noz+1, sizeof(double));
    if(E->refstate.rho == NULL || E->refstate.ala_beta == NULL ||
       E->refstate.ala_beta_supplied == NULL ||
       E->refstate.ala_beta_density == NULL ||
       E->refstate.ala_beta_interval == NULL ||
       E->refstate.beta_ala == NULL) {
        fprintf(stderr, "Unable to allocate reference density/ALA beta storage\n");
        parallel_process_termination();
    }
    E->refstate.has_beta_ala = 0;
    E->refstate.has_Ks = 0;
    E->refstate.has_beta_interval = 0;

    /* reference profile of gravity */
    E->refstate.gravity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of coefficient of thermal expansion */
    E->refstate.thermal_expansivity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of heat capacity */
    E->refstate.heat_capacity = (double *) malloc((noz+1)*sizeof(double));

    // DJB EBA
    /* dissipation scaling */
    E->refstate.dis = (double *) malloc((noz+1)*sizeof(double));

    /* Column 3 is the fixed thermodynamic reference temperature. */
    E->refstate.temperature = (double *) malloc((noz+1)*sizeof(double));
    /* Tref is a semantic alias. temperature remains the storage owner. */
    E->refstate.Tref = E->refstate.temperature;
    E->refstate.gamma_eff = (double *) malloc((noz+1)*sizeof(double));
    E->refstate.Ks = (double *) calloc(noz+1, sizeof(double));
    if(E->refstate.gravity == NULL ||
       E->refstate.thermal_expansivity == NULL ||
       E->refstate.heat_capacity == NULL ||
       E->refstate.dis == NULL ||
       E->refstate.temperature == NULL ||
       E->refstate.gamma_eff == NULL ||
       E->refstate.Ks == NULL) {
        fprintf(stderr, "Unable to allocate reference thermodynamic storage\n");
        parallel_process_termination();
    }
    E->refstate.has_temperature = 0;
    E->refstate.temperature_cmb = 0.0;
    E->refstate.temperature_surface = 0.0;

    /* reference profile of temperature */
    /*E->refstate.Tadi = (double *) malloc((noz+1)*sizeof(double));*/

}


void mat_prop_free(struct All_variables *E)
{
    free(E->refstate.ala_beta);
    E->refstate.ala_beta = NULL;
    free(E->refstate.ala_beta_supplied);
    E->refstate.ala_beta_supplied = NULL;
    free(E->refstate.ala_beta_density);
    E->refstate.ala_beta_density = NULL;
    free(E->refstate.ala_beta_interval);
    E->refstate.ala_beta_interval = NULL;
    free(E->refstate.beta_ala);
    E->refstate.beta_ala = NULL;
    free(E->refstate.Ks);
    E->refstate.Ks = NULL;
}


void reference_state(struct All_variables *E)
{
    int i;

    /* All refstate variables (except Tadi) must be 1 at the surface.
     * Otherwise, the scaling of eqns in the code might not be correct. */

    /* select the choice of reference state */
    switch(E->refstate.choice) {
    case 0:
        /* read from a file */
        read_refstate(E);
        if(strcmp(E->control.ala_beta_element_source,"interval") == 0)
            read_ala_beta_intervals(E);
        break;
    case 1:
        /* Adams-Williamson EoS */
        adams_williamson_eos(E);
        break;
    default:
        if (E->parallel.me) {
            fprintf(stderr, "Unknown option for reference state\n");
            fprintf(E->fp, "Unknown option for reference state\n");
            fflush(E->fp);
        }
        parallel_process_termination();
    }

    initialize_ala_beta(E);

    if(E->control.ala_pressure_buoyancy) {
        if(E->parallel.me == 0) {
            fprintf(stderr, "nz  radius   depth    rho    T_init_background    beta    gamma_eff    layer\n");
        }
        if(E->parallel.me < E->parallel.nprocz)
            for(i=1; i<=E->lmesh.noz; i++) {
                fprintf(stderr, "%d %f %f %e %11f %e %e %5i\n",
                        i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                        E->refstate.rho[i], E->refstate.temperature[i],
                        E->refstate.beta_ala[i], E->refstate.gamma_eff[i],
                        layers_r(E,E->sx[1][3][i]));
            }
    }
    else {
        if(E->parallel.me == 0) {
            fprintf(stderr, "nz  radius   depth    rho    dis     layer\n");  // DJB EBA
        }
        if(E->parallel.me < E->parallel.nprocz)
            for(i=1; i<=E->lmesh.noz; i++) {
            fprintf(stderr, "%d %f %f %e %11f %5i\n",
                    i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                    E->refstate.rho[i],E->refstate.dis[i],layers_r(E,E->sx[1][3][i]));  // DJB EBA
            }
    }

    return;
}


/* Select one authoritative strict-ALA element beta.  Continuity, pressure
 * buoyancy, augmented-Lagrangian terms, multigrid restriction and diagnostics
 * all consume ala_beta.  File-backed interval beta is the strict,
 * phase-inclusive density representation; legacy paths retain nodal endpoint
 * averaging or the density logarithmic secant.
 *
 * Validation uses relative mismatch |beta-beta_rho|/max(|beta_rho|,1e-12).
 * RMS > 2% or maximum > 5% warns; RMS > 10% or maximum > 25% fails.  These
 * thresholds also make missing R0 or kilometre-to-metre conversions fatal.
 */
static void initialize_ala_beta(struct All_variables *E)
{
    int nz, global_nz, local_bad, global_bad, local_count, global_count;
    double r0, r1, rho0, rho1, beta, beta_rho, beta_thermo;
    double mismatch, relative, thermo_mismatch, thermo_relative;
    double local_sum, global_sum, local_sq, global_sq;
    double local_rel_sq, global_rel_sq;
    double local_thermo_sum, global_thermo_sum;
    double local_thermo_sq, global_thermo_sq;
    double local_thermo_rel_sq, global_thermo_rel_sq;
    double local_signed_at_max, global_signed_at_max;
    double local_depth_at_max, global_depth_at_max;
    double rms, relative_rms, thermo_rms, thermo_relative_rms;
    struct { double value; int index; } local_max, global_max;
    struct { double value; int index; } local_thermo_max, global_thermo_max;
    const double beta_floor = 1.0e-12;

    for(nz=1; nz<=E->lmesh.noz; nz++) {
        if(E->refstate.rho[nz] <= 0.0 || !isfinite(E->refstate.rho[nz])) {
            fprintf(stderr,
                    "Invalid smooth reference density: rank=%d global_nz=%d rho=%e\n",
                    E->parallel.me, nz + E->lmesh.nzs - 1,
                    E->refstate.rho[nz]);
            parallel_process_termination();
        }
    }

    local_sum = local_sq = local_rel_sq = 0.0;
    local_thermo_sum = local_thermo_sq = local_thermo_rel_sq = 0.0;
    local_max.value = local_thermo_max.value = -1.0;
    local_max.index = local_thermo_max.index = -1;
    local_bad = 0;
    local_count = E->lmesh.elz;

    for(nz=1; nz<=E->lmesh.elz; nz++) {
        r0 = E->sx[1][3][nz];
        r1 = E->sx[1][3][nz+1];
        rho0 = E->refstate.rho[nz];
        rho1 = E->refstate.rho[nz+1];
        if(!isfinite(r0) || !isfinite(r1) || r1 <= r0) {
            fprintf(stderr,
                    "Invalid outward radial interval: rank=%d local_nz=%d r0=%e r1=%e\n",
                    E->parallel.me, nz, r0, r1);
            parallel_process_termination();
        }

        beta_rho = -(log(rho1) - log(rho0)) / (r1 - r0);
        E->refstate.ala_beta_supplied[nz] =
            0.5 * (E->refstate.beta_ala[nz] +
                   E->refstate.beta_ala[nz+1]);
        E->refstate.ala_beta_density[nz] = beta_rho;
        if(E->refstate.choice != 0)
            E->refstate.ala_beta_supplied[nz] = beta_rho;
        if(E->control.ala_pressure_buoyancy && E->refstate.choice == 0 &&
           strcmp(E->control.ala_beta_element_source,"supplied_average") == 0)
            beta = E->refstate.ala_beta_supplied[nz];
        else if(E->control.ala_pressure_buoyancy &&
                E->refstate.choice == 0 &&
                strcmp(E->control.ala_beta_element_source,"interval") == 0)
            beta = E->refstate.ala_beta_interval[nz];
        else
            beta = E->refstate.ala_beta_density[nz];
        if(!isfinite(E->refstate.ala_beta_supplied[nz]) ||
           E->refstate.ala_beta_supplied[nz] <= 0.0 ||
           !isfinite(E->refstate.ala_beta_density[nz]) ||
           E->refstate.ala_beta_density[nz] <= 0.0 ||
           !isfinite(beta) || beta <= 0.0) {
            fprintf(stderr,
                    "Invalid ALA beta: rank=%d local_nz=%d "
                    "supplied=%e density=%e selected=%e; "
                    "strict ALA requires positive beta*=-d(ln rho)/dr*\n",
                    E->parallel.me, nz,
                    E->refstate.ala_beta_supplied[nz],
                    E->refstate.ala_beta_density[nz],beta);
            parallel_process_termination();
        }
        E->refstate.ala_beta[nz] = beta;

        mismatch = beta - beta_rho;
        relative = fabs(mismatch) / max(fabs(beta_rho), beta_floor);
        beta_thermo = E->control.disptn_number
            * 0.5 * (E->refstate.thermal_expansivity[nz] +
                     E->refstate.thermal_expansivity[nz+1])
            * 0.5 * (E->refstate.gravity[nz] + E->refstate.gravity[nz+1])
            / (0.5 * (E->refstate.heat_capacity[nz] +
                      E->refstate.heat_capacity[nz+1])
               * 0.5 * (E->refstate.gamma_eff[nz] +
                        E->refstate.gamma_eff[nz+1]));
        thermo_mismatch = beta - beta_thermo;
        thermo_relative = fabs(thermo_mismatch) /
                          max(fabs(beta_thermo), beta_floor);
        global_nz = nz + E->lmesh.nzs - 1;

        local_sum += mismatch;
        local_sq += mismatch * mismatch;
        local_rel_sq += relative * relative;
        local_thermo_sum += thermo_mismatch;
        local_thermo_sq += thermo_mismatch * thermo_mismatch;
        local_thermo_rel_sq += thermo_relative * thermo_relative;
        if(relative > local_max.value) {
            local_max.value = relative;
            local_max.index = global_nz;
        }
        if(thermo_relative > local_thermo_max.value) {
            local_thermo_max.value = thermo_relative;
            local_thermo_max.index = global_nz;
        }
        if(relative > 0.25 || thermo_relative > 0.25)
            local_bad = 1;
    }

    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_sq, &global_sq, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_rel_sq, &global_rel_sq, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_thermo_sum, &global_thermo_sum, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_thermo_sq, &global_thermo_sq, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_thermo_rel_sq, &global_thermo_rel_sq, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, E->parallel.world);
    MPI_Allreduce(&local_thermo_max, &global_thermo_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, E->parallel.world);
    MPI_Allreduce(&local_bad, &global_bad, 1, MPI_INT, MPI_MAX, E->parallel.world);
    rms = sqrt(global_sq / global_count);
    relative_rms = sqrt(global_rel_sq / global_count);
    thermo_rms = sqrt(global_thermo_sq / global_count);
    thermo_relative_rms = sqrt(global_thermo_rel_sq / global_count);

    local_signed_at_max = local_depth_at_max = 0.0;
    for(nz=1; nz<=E->lmesh.elz; nz++) {
        global_nz = nz + E->lmesh.nzs - 1;
        if(global_nz == global_max.index) {
            beta_rho = -(log(E->refstate.rho[nz+1]) - log(E->refstate.rho[nz])) /
                       (E->sx[1][3][nz+1] - E->sx[1][3][nz]);
            local_signed_at_max = E->refstate.ala_beta[nz] - beta_rho;
            local_depth_at_max = 1.0 - 0.5 * (E->sx[1][3][nz] + E->sx[1][3][nz+1]);
        }
    }
    MPI_Allreduce(&local_signed_at_max, &global_signed_at_max, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    MPI_Allreduce(&local_depth_at_max, &global_depth_at_max, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);
    /* Horizontal ranks duplicate the same radial reference-state slice. */
    global_signed_at_max /= E->parallel.nprocxy;
    global_depth_at_max /= E->parallel.nprocxy;

    if(E->parallel.me == 0 && E->control.ala_pressure_buoyancy) {
        fprintf(stderr,
                "Strict ALA beta validation source=%s "
                "(beta* uses radius normalized by R0)\n"
                "  density closure: mean_signed=%e RMS=%e relative_RMS=%e max_relative=%e "
                "signed_at_max=%e radial_element=%d depth_over_R0=%e\n"
                "  thermodynamic closure: mean_signed=%e RMS=%e "
                "relative_RMS=%e max_relative=%e radial_element=%d\n",
                E->control.ala_beta_element_source,
                global_sum/global_count, rms, relative_rms, global_max.value,
                global_signed_at_max, global_max.index, global_depth_at_max,
                global_thermo_sum/global_count, thermo_rms, thermo_relative_rms,
                global_thermo_max.value, global_thermo_max.index);
        if(relative_rms > 0.02 || thermo_relative_rms > 0.02 ||
           global_max.value > 0.05 || global_thermo_max.value > 0.05)
            fprintf(stderr, "  WARNING: strict ALA closure exceeds 2%% RMS or 5%% maximum tolerance\n");
    }
    if(E->control.ala_pressure_buoyancy &&
       (global_bad || relative_rms > 0.10 || thermo_relative_rms > 0.10)) {
        if(E->refstate.has_Ks) {
            if(E->parallel.me == 0)
                fprintf(stderr,
                        "WARNING: phase-inclusive 8-column strict reference "
                        "state exceeds legacy beta-closure tolerance; "
                        "diagnostic only, supplied beta is unchanged\n");
        }
        else {
            fprintf(stderr,
                    "Strict ALA reference state fails beta closure tolerance\n");
            parallel_process_termination();
        }
    }
}


double conductivity_depth_factor(struct All_variables *E,
                                 double physical_depth_m)
{
    double depth, d;

    if(!isfinite(physical_depth_m)) {
        fprintf(stderr, "Non-finite physical depth in conductivity_depth_factor\n");
        parallel_process_termination();
    }

    depth = max(0.0, physical_depth_m);
    d = depth / (E->control.kd_mantle_thickness_km * 1.0e3);
    if(depth < E->control.kd_transition_depth_km * 1.0e3)
        return E->control.kd_upper_prefactor
               * (1.0 + E->control.kd_upper_linear*d
                  + E->control.kd_upper_quadratic*d*d)
               / E->data.ks;
    return E->control.kd_lower_prefactor
           * (1.0 + E->control.kd_lower_linear*d
              + E->control.kd_lower_quadratic*d*d)
           / E->data.ks;
}


double conductivity_temperature_factor(struct All_variables *E,
                                       double nondimensional_temperature)
{
    double dimensional_temperature, surface_temperature;

    if(E->control.kT_exponent == 0.0)
        return 1.0;

    surface_temperature = E->data.Ttop;
    dimensional_temperature = surface_temperature
        + nondimensional_temperature * E->data.ref_temperature;
    /* Do not let numerical undershoot alter the evolved T*.  The conductivity
       lookup alone is bounded at the explicit dimensional top temperature. */
    if(dimensional_temperature < surface_temperature)
        dimensional_temperature = surface_temperature;

    return pow(surface_temperature / dimensional_temperature,
               (double)E->control.kT_exponent);
}


static double conductivity_composition_factor(struct All_variables *E,
                                              double primordial_fraction)
{
    double cprim;

    if(!isfinite(primordial_fraction)) {
        fprintf(stderr, "Non-finite primordial fraction in conductivity model\n");
        parallel_process_termination();
    }
    cprim = min(1.0, max(0.0, primordial_fraction));
    return 1.0 + (E->control.kC_ratio - 1.0) * cprim;
}


double conductivity_element_composition_factor(struct All_variables *E,
                                               int cap, int element)
{
    /* The ratio composition scheme maps component i to tracer flavor i+1.
       The highest-numbered flavor is primordial in both supported setups. */
    if(!E->composition.on || E->composition.ncomp < 1)
        return 1.0;
    return conductivity_composition_factor(
        E, E->composition.comp_el[cap][E->composition.ncomp-1][element]);
}


double conductivity_element_prefactor(struct All_variables *E,
                                      int cap, int element,
                                      double reference_conductivity)
{
    double depth_m, kd, kC;

    depth_m = (1.0 - 0.5
        * (E->sx[cap][3][E->ien[cap][element].node[1]]
           + E->sx[cap][3][E->ien[cap][element].node[5]]))
        * E->data.radius_km * 1.0e3;
    kd = conductivity_depth_factor(E, depth_m);
    kC = conductivity_element_composition_factor(E, cap, element);
    return reference_conductivity * kd * kC;
}


double nodal_conductivity_diagnostic(struct All_variables *E, int cap, int node,
                                    int component)
{
    int iz;
    double depth_m, kd, kT, kC, k_total, rho, cp;

    iz = ((node - 1) % E->lmesh.noz) + 1;
    depth_m = (1.0 - E->sx[cap][3][node]) * E->data.radius_km * 1.0e3;
    kd = conductivity_depth_factor(E, depth_m);
    kT = conductivity_temperature_factor(E, E->T[cap][node]);
    if(E->composition.on && E->composition.ncomp > 0)
        kC = conductivity_composition_factor(
            E, E->composition.comp_node[cap][E->composition.ncomp-1][node]);
    else
        kC = 1.0;
    k_total = kd * kT * kC;
    rho = E->refstate.rho[iz];
    cp = E->refstate.heat_capacity[iz];

    switch(component) {
    case CONDUCTIVITY_KD: return kd;
    case CONDUCTIVITY_KT: return kT;
    case CONDUCTIVITY_KC: return kC;
    case CONDUCTIVITY_K_TOTAL:
        return E->control.reference_conductivity * k_total;
    case CONDUCTIVITY_KAPPA_EFF:
        return E->control.reference_conductivity * k_total / (rho * cp);
    case CONDUCTIVITY_RHO_REF: return rho;
    case CONDUCTIVITY_CP: return cp;
    default:
        fprintf(stderr, "Unknown conductivity diagnostic component %d\n",
                component);
        parallel_process_termination();
    }
    return 0.0;
}


double nodal_thermal_conductivity(struct All_variables *E, int cap, int node)
{
    return nodal_conductivity_diagnostic(E, cap, node,
                                         CONDUCTIVITY_K_TOTAL);
}


static void read_refstate(struct All_variables *E)
{
    FILE *fp;
    int i, columns, expected_columns, cmb_columns, background_rows;
    char buffer[255], cmb_buffer[255], trailing;
    double values[9], cmb_values[9];
    double last_temperature;

    fp = fopen(E->refstate.filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open reference state file: %s\n",
                E->refstate.filename);
        parallel_process_termination();
    }

    /* Every radial MPI rank needs the same reference CMB temperature
       when constructing the bottom TBL. Read the first global row explicitly,
       then rewind before the existing local radial slice read. */
    if(!read_refstate_data_line(fp, cmb_buffer, 255)) {
        fprintf(stderr, "Reference state file '%s' is empty\n",
                E->refstate.filename);
        parallel_process_termination();
    }
    if(E->control.ala_pressure_buoyancy)
        cmb_columns = sscanf(cmb_buffer,
                             "%lf %lf %lf %lf %lf %lf %lf %lf %c",
                             &cmb_values[0], &cmb_values[1], &cmb_values[2],
                             &cmb_values[3], &cmb_values[4], &cmb_values[5],
                             &cmb_values[6], &cmb_values[7], &trailing);
    else
        cmb_columns = sscanf(cmb_buffer,
                             "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                             &cmb_values[0], &cmb_values[1], &cmb_values[2],
                             &cmb_values[3], &cmb_values[4], &cmb_values[5],
                             &cmb_values[6], &cmb_values[7], &cmb_values[8]);
    if(E->control.ala_pressure_buoyancy &&
       cmb_columns != 7 && cmb_columns != 8) {
        fprintf(stderr,
                "Reference state file '%s', global radial row 1: "
                "strict ALA requires 7 or 8 numeric columns "
                "(rho g Tref alpha Cp beta Gamma_eff [Ks]), found %d\n",
                E->refstate.filename, cmb_columns);
        parallel_process_termination();
    }
    if(!E->control.ala_pressure_buoyancy &&
       cmb_columns != 7 && cmb_columns != 8 && cmb_columns != 9) {
        fprintf(stderr,
                "Legacy reference state file '%s', global radial row 1: "
                "expected 7, 8, or 9 numeric columns, found %d\n",
                E->refstate.filename, cmb_columns);
        parallel_process_termination();
    }
    if(E->control.ala_pressure_buoyancy) {
        /* Column 3 is the thermodynamic Tref on every node, including the
         * endpoints. Boundary temperatures remain part of the total-
         * temperature initialization and are not encoded in Tref. */
        rewind(fp);
        background_rows = 0;
        last_temperature = 0.0;
        while(read_refstate_data_line(fp, buffer, 255)) {
            columns = sscanf(buffer,
                             "%lf %lf %lf %lf %lf %lf %lf %lf %c",
                             &values[0], &values[1], &values[2], &values[3],
                             &values[4], &values[5], &values[6], &values[7],
                             &trailing);
            if(columns != 7 && columns != 8) {
                fprintf(stderr,
                        "Reference state file '%s', global radial row %d: "
                        "strict ALA requires 7 or 8 numeric columns\n",
                        E->refstate.filename, background_rows+1);
                parallel_process_termination();
            }
            if(background_rows == 0)
                E->refstate.temperature_cmb = values[2];
            last_temperature = values[2];
            background_rows++;
        }
        if(background_rows != E->mesh.noz) {
            fprintf(stderr,
                    "Reference state file '%s' has %d radial rows; expected %d\n",
                    E->refstate.filename, background_rows, E->mesh.noz);
            parallel_process_termination();
        }
        E->refstate.temperature_surface = last_temperature;
    }
    else if(cmb_columns >= 8)
        E->refstate.temperature_cmb = cmb_values[6];
    rewind(fp);

    /* skip these lines, which belong to other processors */
    for(i=1; i<E->lmesh.nzs; i++) {
        if(!read_refstate_data_line(fp, buffer, 255)) {
            fprintf(stderr, "Reference state file '%s' has too few radial rows\n",
                    E->refstate.filename);
            parallel_process_termination();
        }
    }

    expected_columns = 0;
    for(i=1; i<=E->lmesh.noz; i++) {
        if(!read_refstate_data_line(fp, buffer, 255)) {
            fprintf(stderr, "Reference state file '%s' has too few radial rows\n",
                    E->refstate.filename);
            parallel_process_termination();
        }

        if(E->control.ala_pressure_buoyancy)
            columns = sscanf(buffer,
                             "%lf %lf %lf %lf %lf %lf %lf %lf %c",
                             &values[0], &values[1], &values[2], &values[3],
                             &values[4], &values[5], &values[6], &values[7],
                             &trailing);
        else
            columns = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                             &values[0], &values[1], &values[2], &values[3],
                             &values[4], &values[5], &values[6], &values[7],
                             &values[8]);
        if(E->control.ala_pressure_buoyancy &&
           columns != 7 && columns != 8) {
            fprintf(stderr,
                    "Reference state file '%s', global radial row %d: "
                    "strict ALA requires 7 or 8 numeric columns "
                    "(rho g Tref alpha Cp beta Gamma_eff [Ks]), found %d\n",
                    E->refstate.filename, i+E->lmesh.nzs-1, columns);
            parallel_process_termination();
        }
        if(!E->control.ala_pressure_buoyancy &&
           columns != 7 && columns != 8 && columns != 9) {
            fprintf(stderr,
                    "Legacy reference state file '%s', global radial row %d: "
                    "expected 7, 8, or 9 numeric columns, found %d\n",
                    E->refstate.filename, i+E->lmesh.nzs-1, columns);
            parallel_process_termination();
        }
        if(expected_columns == 0)
            expected_columns = columns;
        else if(columns != expected_columns) {
            fprintf(stderr,
                    "Reference state file '%s' mixes %d- and %d-column rows\n",
                    E->refstate.filename, expected_columns, columns);
            parallel_process_termination();
        }

        if(E->control.ala_pressure_buoyancy) {
            /* Strict schema:
             * rho g Tref alpha Cp beta Gamma_eff [Ks]. */
            E->refstate.rho[i] = values[0];
            E->refstate.gravity[i] = values[1];
            E->refstate.temperature[i] = values[2];
            E->refstate.thermal_expansivity[i] = values[3];
            E->refstate.heat_capacity[i] = values[4];
            E->refstate.beta_ala[i] = values[5];
            E->refstate.gamma_eff[i] = values[6];
            E->refstate.Ks[i] = columns == 8 ? values[7] : 0.0;
            E->refstate.dis[i] = 1.0; /* inert legacy storage, not strict input */
            if(E->refstate.rho[i] <= 0.0 ||
               E->refstate.gravity[i] <= 0.0 ||
               E->refstate.thermal_expansivity[i] <= 0.0 ||
               E->refstate.heat_capacity[i] <= 0.0 ||
               E->refstate.beta_ala[i] <= 0.0 ||
               E->refstate.gamma_eff[i] <= 0.0 ||
               !isfinite(E->refstate.rho[i]) ||
               !isfinite(E->refstate.gravity[i]) ||
               !isfinite(E->refstate.temperature[i]) ||
               !isfinite(E->refstate.thermal_expansivity[i]) ||
               !isfinite(E->refstate.heat_capacity[i]) ||
               !isfinite(E->refstate.beta_ala[i]) ||
               !isfinite(E->refstate.gamma_eff[i]) ||
               (columns == 8 &&
                (!isfinite(E->refstate.Ks[i]) ||
                 E->refstate.Ks[i] <= 0.0))) {
                fprintf(stderr,
                        "Invalid strict ALA row %d: all fields must be finite "
                        "and rho/g/alpha/Cp/beta/Gamma_eff plus optional Ks "
                        "positive\n",
                        i+E->lmesh.nzs-1);
                parallel_process_termination();
            }
        }
        else {
            /* Historical named formulations retain their dis-containing
             * rho g alpha Cp dis k [Tref Gamma_eff [beta]] interpretation. */
            E->refstate.rho[i] = values[0];
            E->refstate.gravity[i] = values[1];
            E->refstate.thermal_expansivity[i] = values[2];
            E->refstate.heat_capacity[i] = values[3];
            E->refstate.dis[i] = values[4];
            E->refstate.temperature[i] = columns >= 8 ? values[6] : 0.0;
            E->refstate.gamma_eff[i] = columns >= 8 ? values[7] : 1.0;
            E->refstate.beta_ala[i] = columns == 9 ? values[8] : 0.0;
            E->refstate.Ks[i] = 0.0;
        }

        /**** debug ****
        fprintf(stderr, "%d %f %f %f %f\n",
                i,
                E->refstate.rho[i],
                E->refstate.gravity[i],
                E->refstate.thermal_expansivity[i],
                E->refstate.heat_capacity[i]);
        /* end of debug */
    }

    E->refstate.has_Ks = E->control.ala_pressure_buoyancy &&
                         expected_columns == 8;
    E->refstate.has_temperature = E->control.ala_pressure_buoyancy ||
                                  expected_columns >= 8;
    E->refstate.has_beta_ala = E->control.ala_pressure_buoyancy ||
                               expected_columns == 9;
    if(!E->refstate.has_temperature && E->control.lith_age &&
       E->convection.tic_method != -1) {
        fprintf(stderr,
                "Initial-background geotherm construction requires an 8-column "
                "reference state file; '%s' uses the legacy 7-column format\n",
                E->refstate.filename);
        parallel_process_termination();
    }

    if(E->parallel.me == 0 && E->control.ala_pressure_buoyancy) {
        fprintf(stderr,
                "Read strict ALA reference state '%s': "
                "rho g Tref alpha Cp beta Gamma_eff%s; "
                "Tref endpoints CMB=%e surface=%e\n",
                E->refstate.filename,
                expected_columns == 8 ? " Ks" : "",
                E->refstate.temperature_cmb,
                E->refstate.temperature_surface);
    }
    else if(E->parallel.me == 0 && E->refstate.has_temperature) {
        fprintf(stderr,
                "Read %d-column legacy dis-containing reference state '%s'\n",
                expected_columns, E->refstate.filename);
    }

    fclose(fp);
    return;
}


static void read_ala_beta_intervals(struct All_variables *E)
{
    FILE *fp;
    int columns, element_index, expected_index, local_index;
    int local_first, local_last;
    char buffer[255], trailing;
    double r_inner, r_outer, beta_interval;
    double mesh_inner, mesh_outer;
    const double radius_tolerance = 1.0e-10;

    if(!E->refstate.has_Ks) {
        fprintf(stderr,
                "ala_beta_element_source=interval requires an 8-column "
                "strict reference state\n");
        parallel_process_termination();
    }

    fp = fopen(E->refstate.beta_interval_filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "Cannot open strict ALA interval beta file: %s\n",
                E->refstate.beta_interval_filename);
        parallel_process_termination();
    }

    expected_index = 1;
    local_first = E->lmesh.nzs;
    local_last = local_first + E->lmesh.elz - 1;
    while(read_refstate_data_line(fp, buffer, 255)) {
        columns = sscanf(buffer, "%d %lf %lf %lf %c",
                         &element_index, &r_inner, &r_outer,
                         &beta_interval, &trailing);
        if(columns != 4 || element_index != expected_index) {
            fprintf(stderr,
                    "Interval beta file '%s', row %d: expected "
                    "'element_index r_inner r_outer beta_interval' with "
                    "sequential index %d\n",
                    E->refstate.beta_interval_filename, expected_index,
                    expected_index);
            parallel_process_termination();
        }
        if(!isfinite(r_inner) || !isfinite(r_outer) ||
           !isfinite(beta_interval) || r_outer <= r_inner ||
           beta_interval <= 0.0) {
            fprintf(stderr,
                    "Invalid interval beta row %d in '%s'\n",
                    element_index, E->refstate.beta_interval_filename);
            parallel_process_termination();
        }

        if(element_index >= local_first && element_index <= local_last) {
            local_index = element_index - local_first + 1;
            mesh_inner = E->sx[1][3][local_index];
            mesh_outer = E->sx[1][3][local_index+1];
            if(fabs(r_inner-mesh_inner) > radius_tolerance ||
               fabs(r_outer-mesh_outer) > radius_tolerance) {
                fprintf(stderr,
                        "Interval beta mesh mismatch: rank=%d global_element=%d "
                        "file=(%.16e,%.16e) mesh=(%.16e,%.16e)\n",
                        E->parallel.me, element_index, r_inner, r_outer,
                        mesh_inner, mesh_outer);
                parallel_process_termination();
            }
            E->refstate.ala_beta_interval[local_index] = beta_interval;
        }
        expected_index++;
    }
    fclose(fp);

    if(expected_index-1 != E->mesh.elz) {
        fprintf(stderr,
                "Interval beta file '%s' has %d rows; expected %d radial "
                "elements\n",
                E->refstate.beta_interval_filename, expected_index-1,
                E->mesh.elz);
        parallel_process_termination();
    }
    for(local_index=1; local_index<=E->lmesh.elz; local_index++) {
        if(!isfinite(E->refstate.ala_beta_interval[local_index]) ||
           E->refstate.ala_beta_interval[local_index] <= 0.0) {
            fprintf(stderr,
                    "Missing interval beta: rank=%d local_element=%d\n",
                    E->parallel.me, local_index);
            parallel_process_termination();
        }
    }
    E->refstate.has_beta_interval = 1;

    if(E->parallel.me == 0)
        fprintf(stderr,
                "Read strict ALA interval beta '%s': %d radial elements\n",
                E->refstate.beta_interval_filename, E->mesh.elz);
}


static int read_refstate_data_line(FILE *fp, char *buffer, int length)
{
    char *cursor;

    while(fgets(buffer, length, fp) != NULL) {
        cursor = buffer;
        while(isspace((unsigned char)*cursor))
            cursor++;
        if(*cursor != '\0' && *cursor != '#')
            return 1;
    }
    return 0;
}


static void adams_williamson_eos(struct All_variables *E)
{
    int i;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    for(i=1; i<=E->lmesh.noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.rho[i] = exp(beta*z);
	E->refstate.gravity[i] = 1;
	E->refstate.thermal_expansivity[i] = 1;
	E->refstate.heat_capacity[i] = 1;
	E->refstate.dis[i] = 1; // DJB EBA
	E->refstate.temperature[i] = 0.0;
	E->refstate.gamma_eff[i] = 1.0;
	E->refstate.Ks[i] = 0.0;
	/*E->refstate.Tadi[i] = (E->control.adiabaticT0 + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;*/
    }

    E->refstate.temperature_cmb = 0.0;
    E->refstate.temperature_surface = 0.0;
    E->refstate.has_beta_ala = 0;
    E->refstate.has_Ks = 0;
    for(i=1; i<=E->lmesh.noz; i++)
        E->refstate.beta_ala[i] = 0.0;

    return;
}
