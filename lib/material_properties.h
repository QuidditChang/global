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


void mat_prop_allocate(struct All_variables *E);
void mat_prop_free(struct All_variables *E);
void reference_state(struct All_variables *E);
double conductivity_depth_factor(struct All_variables *E,
                                 double physical_depth_m);
double conductivity_temperature_factor(struct All_variables *E,
                                       double nondimensional_temperature);
double conductivity_element_composition_factor(struct All_variables *E,
                                               int cap, int element);
double conductivity_element_prefactor(struct All_variables *E,
                                      int cap, int element,
                                      double reference_conductivity);
double nodal_thermal_conductivity(struct All_variables *E, int cap, int node);
double nodal_conductivity_diagnostic(struct All_variables *E, int cap, int node,
                                    int component);

#define CONDUCTIVITY_KD 0
#define CONDUCTIVITY_KT 1
#define CONDUCTIVITY_KC 2
#define CONDUCTIVITY_K_TOTAL 3
#define CONDUCTIVITY_KAPPA_EFF 4
#define CONDUCTIVITY_RHO_REF 5
#define CONDUCTIVITY_CP 6


void density(struct All_variables *E, double *rho);
