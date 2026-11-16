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
void phase_change_input(struct All_variables *E);
void phase_change_allocate(struct All_variables *E);

void phase_change_apply(struct All_variables *E, double **buoy);
int phase_change_count(const struct All_variables *E);
double phase_change_entropy(const struct All_variables *E, int phase_index);
int phase_change_thermal_active(const struct All_variables *E, int phase_index);
void phase_change_set_pPv(struct All_variables *E, float depth, float transT,
                          float clapeyron, float width, float density_jump,
                          float entropy_jump);
/* fraction is the already evaluated high-pressure phase fraction, not T. */
double phase_change_density_fraction(const struct All_variables *E,
                                      int phase_index, int cap, int node,
                                      double fraction);
double pPv_fraction_at_state(const struct All_variables *E, double depth,
                             double temperature, double rho, double gravity);
double pPv_fraction_at_node(const struct All_variables *E, int cap, int node,
                            double temperature);
double pPv_fraction_at_gp(const struct All_variables *E, int cap, int element,
                          int gp, double temperature);
double pPv_conductivity_factor(const struct All_variables *E, double fraction);
double pPv_viscosity_factor(const struct All_variables *E, double fraction);
void phase_change_state(const struct Phase_transition *phase,
                        double depth, double temperature,
                        double rho, double gravity,
                        double d_rho_g_dr,
                        double *q, double *fraction,
                        double *dX_dT, double *dX_dr_pressure);
float phase_change_reference_fraction(const struct All_variables *E,
                                      int phase_index, int cap, int node);
float phase_change_fraction_at_temperature(const struct All_variables *E,
                                           int phase_index, int cap, int node,
                                           double temperature);
