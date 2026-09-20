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
#if !defined(CitcomS_drive_solvers_h)
#define CitcomS_drive_solvers_h

#ifdef __cplusplus
extern "C" {
#endif

/* Cached, pre-rigid-rotation EBA Stokes snapshot. Powers use Di/Atemp;
 * shell integrals have the same units as totals, ordered surface to CMB. */
enum Eba_power_term {
    EBA_PPLATE, EBA_POTHER, EBA_WTRACTION, EBA_WBODY, EBA_DOPERATOR,
    EBA_QVISC, EBA_WTHERMAL, EBA_WCHEMICAL, EBA_W410, EBA_W520,
    EBA_W660, EBA_WPHASE, EBA_WPRESSURE, EBA_RMECHANICAL,
    EBA_ROPERATOR, EBA_RBODY_SPLIT, EBA_RHEATING_OPERATOR,
    EBA_QVISC_CAPPED, EBA_QVISC_USED, EBA_QVISC_REMOVED,
    EBA_QVISC_POTENTIAL, EBA_QVISC_LIMITED_VOLUME, EBA_POWER_COUNT
};
struct Eba_power_snapshot {
    int step, depth_count, qvis_mode;
    double elapsed_time, scale;
    double total[EBA_POWER_COUNT];
    double *shell_integrals;
};
extern const char *const eba_power_names[EBA_POWER_COUNT];
const struct Eba_power_snapshot *eba_power_snapshot(struct All_variables *E);

void general_stokes_solver(struct All_variables*);
void general_stokes_solver_setup(struct All_variables*);

#ifdef __cplusplus
}
#endif

#endif
