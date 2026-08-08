/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//<LicenseText>
//
// CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
// Copyright (C) 2002-2005, California Institute of Technology.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//</LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#include <Python.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global_defs.h"
#include "parallel_related.h"
#include "setProperties.h"


/* See PEP 353. */
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif


/*==============================================================*/
/* functions and macros which fetch properties from 'inventory' */


FILE *get_output_stream(PyObject *out, struct All_variables*E);
#define PUTS(s) if (fp) fprintf(fp, s)

int _getStringProperty(PyObject* properties, char* attribute,
                       char* value, size_t valueSize, FILE* fp);
#define getStringProperty(p, a, v, o) if (-1 == _getStringProperty(p, a, v, sizeof(v), o)) return NULL

int _getIntProperty(PyObject* properties, char* attribute, int *value, FILE* fp);
#define getIntProperty(p, a, v, o) if (-1 == _getIntProperty(p, a, &(v), o)) return NULL

int _getFloatProperty(PyObject* properties, char* attribute, float *value, FILE* fp);
#define getFloatProperty(p, a, v, o) if (-1 == _getFloatProperty(p, a, &(v), o)) return NULL

int _getDoubleProperty(PyObject* properties, char* attribute, double *value, FILE* fp);
#define getDoubleProperty(p, a, v, o) if (-1 == _getDoubleProperty(p, a, &(v), o)) return NULL

int _getIntVectorProperty(PyObject* properties, char* attribute,
                          int* vector, int len, FILE* fp);
#define getIntVectorProperty(p, a, v, l, o) if (-1 == _getIntVectorProperty(p, a, v, l, o)) return NULL

int _getFloatVectorProperty(PyObject* properties, char* attribute,
                            float* vector, int len, FILE* fp);
#define getFloatVectorProperty(p, a, v, l, o) if (-1 == _getFloatVectorProperty(p, a, v, l, o)) return NULL

int _getDoubleVectorProperty(PyObject* properties, char* attribute,
                             double* vector, int len, FILE* fp);
#define getDoubleVectorProperty(p, a, v, l, o) if (-1 == _getDoubleVectorProperty(p, a, v, l, o)) return NULL


void myerror(struct All_variables *,char *);
void report(struct All_variables *,char *);

/*==============================================================*/


char pyCitcom_Advection_diffusion_set_properties__doc__[] = "";
char pyCitcom_Advection_diffusion_set_properties__name__[] = "Advection_diffusion_set_properties";

PyObject * pyCitcom_Advection_diffusion_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    float legacy_inputdiffusivity;
    float reference_conductivity;

    if (!PyArg_ParseTuple(args, "OOO:Advection_diffusion_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.tsolver]\n"));

    getIntProperty(properties, "ADV", E->advection.ADVECTION, fp);
    getIntProperty(properties, "filter_temp", E->advection.filter_temperature, fp);
    getIntProperty(properties, "monitor_max_T", E->advection.monitor_max_T, fp);

    getFloatProperty(properties, "finetunedt", E->advection.fine_tune_dt, fp);
    getFloatProperty(properties, "fixed_timestep", E->advection.fixed_timestep, fp);
    getFloatProperty(properties, "adv_gamma", E->advection.gamma, fp);
    getIntProperty(properties, "adv_sub_iterations", E->advection.temp_iterations, fp);

    getFloatProperty(properties, "inputdiffusivity", legacy_inputdiffusivity, fp);
    getFloatProperty(properties, "reference_conductivity",
                     reference_conductivity, fp);
    if(reference_conductivity >= 0.0)
        E->control.requested_reference_conductivity = reference_conductivity;
    else if(reference_conductivity != -1.0) {
        PyErr_SetString(PyExc_ValueError,
                        "reference_conductivity must be non-negative");
        return NULL;
    }
    else if(legacy_inputdiffusivity >= 0.0)
        E->control.requested_reference_conductivity = legacy_inputdiffusivity;
    else if(legacy_inputdiffusivity == -1.0)
        E->control.requested_reference_conductivity = -1.0;
    else {
        PyErr_SetString(PyExc_ValueError,
                        "inputdiffusivity compatibility value must be non-negative");
        return NULL;
    }
    if(E->data.k0 > 0.0 &&
       E->control.requested_reference_conductivity >= 0.0 &&
       fabs(E->control.requested_reference_conductivity
            - E->control.reference_conductivity)
       > 1.0e-6 * fmax(1.0, fabs(E->control.reference_conductivity))) {
        PyErr_SetString(PyExc_ValueError,
                        "reference_conductivity is derived as ks/k0");
        return NULL;
    }
    getFloatProperty(properties, "kT_exponent", E->control.kT_exponent, fp);
    getFloatProperty(properties, "kC_ratio", E->control.kC_ratio, fp);
    getIntProperty(properties, "kC_primordial_flavor",
                   E->control.kC_primordial_flavor, fp);
    getFloatProperty(properties, "kd_mantle_thickness_km",
                     E->control.kd_mantle_thickness_km, fp);
    getFloatProperty(properties, "kd_transition_depth_km",
                     E->control.kd_transition_depth_km, fp);
    getFloatProperty(properties, "kd_upper_prefactor",
                     E->control.kd_upper_prefactor, fp);
    getFloatProperty(properties, "kd_upper_linear",
                     E->control.kd_upper_linear, fp);
    getFloatProperty(properties, "kd_upper_quadratic",
                     E->control.kd_upper_quadratic, fp);
    getFloatProperty(properties, "kd_lower_prefactor",
                     E->control.kd_lower_prefactor, fp);
    getFloatProperty(properties, "kd_lower_linear",
                     E->control.kd_lower_linear, fp);
    getFloatProperty(properties, "kd_lower_quadratic",
                     E->control.kd_lower_quadratic, fp);


    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;

}



char pyCitcom_BC_set_properties__doc__[] = "";
char pyCitcom_BC_set_properties__name__[] = "BC_set_properties";

PyObject * pyCitcom_BC_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "OOO:BC_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.bc]\n"));

    getIntProperty(properties, "side_sbcs", E->control.side_sbcs, fp);
    getIntProperty(properties, "pseudo_free_surf", E->control.pseudo_free_surf, fp);

    getIntProperty(properties, "topvbc", E->mesh.topvbc, fp);
    getFloatProperty(properties, "topvbxval", E->control.VBXtopval, fp);
    getFloatProperty(properties, "topvbyval", E->control.VBYtopval, fp);

    getIntProperty(properties, "botvbc", E->mesh.botvbc, fp);
    getFloatProperty(properties, "botvbxval", E->control.VBXbotval, fp);
    getFloatProperty(properties, "botvbyval", E->control.VBYbotval, fp);

    getIntProperty(properties, "toptbc", E->mesh.toptbc, fp);
    getFloatProperty(properties, "toptbcval", E->control.TBCtopval, fp);

    getIntProperty(properties, "bottbc", E->mesh.bottbc, fp);
    getFloatProperty(properties, "bottbcval", E->control.TBCbotval, fp);

    getIntProperty(properties, "temperature_bound_adj", E->control.temperature_bound_adj, fp);
    getFloatProperty(properties, "depth_bound_adj", E->control.depth_bound_adj, fp);
    getFloatProperty(properties, "width_bound_adj", E->control.width_bound_adj, fp);


    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_Const_set_properties__doc__[] = "";
char pyCitcom_Const_set_properties__name__[] = "Const_set_properties";

PyObject * pyCitcom_Const_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    float radius;

    if (!PyArg_ParseTuple(args, "OOO:Const_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.const]\n"));

    getFloatProperty(properties, "radius", radius, fp);
    getFloatProperty(properties, "density", E->data.density, fp);
    getFloatProperty(properties, "thermdiff", E->data.therm_diff, fp);
    getFloatProperty(properties, "ks", E->data.ks, fp);
    getFloatProperty(properties, "k0", E->data.k0, fp);
    getFloatProperty(properties, "rho0", E->data.rho0, fp);
    getFloatProperty(properties, "Cp0", E->data.Cp0, fp);
    getFloatProperty(properties, "g0", E->data.g0, fp);
    getFloatProperty(properties, "alpha0", E->data.alpha0, fp);
    getFloatProperty(properties, "Ttop", E->data.Ttop, fp);
    getFloatProperty(properties, "Tbottom", E->data.Tbottom, fp);
    getFloatProperty(properties, "deltaT", E->data.ref_temperature, fp);
    getFloatProperty(properties, "gravacc", E->data.grav_acc, fp);
    getFloatProperty(properties, "thermexp", E->data.therm_exp, fp);
    getFloatProperty(properties, "refvisc", E->data.ref_viscosity, fp);
    getFloatProperty(properties, "cp", E->data.Cp, fp);
    getFloatProperty(properties, "density_above", E->data.density_above, fp);
    getFloatProperty(properties, "density_below", E->data.density_below, fp);

    if(E->data.rho0 <= 0.0)
        E->data.rho0 = E->data.density;
    if(E->data.Cp0 <= 0.0)
        E->data.Cp0 = E->data.Cp;
    if(E->data.g0 <= 0.0)
        E->data.g0 = E->data.grav_acc;
    if(E->data.alpha0 <= 0.0)
        E->data.alpha0 = E->data.therm_exp;
    if(E->data.k0 <= 0.0)
        E->data.k0 = E->data.therm_diff * E->data.rho0 * E->data.Cp0;
    if(E->data.ks <= 0.0)
        E->data.ks = E->data.k0;

    if(!isfinite(E->data.ks) || !isfinite(E->data.k0) ||
       !isfinite(E->data.rho0) || !isfinite(E->data.Cp0) ||
       !isfinite(E->data.g0) || !isfinite(E->data.alpha0) ||
       E->data.ks <= 0.0 || E->data.k0 <= 0.0 ||
       E->data.rho0 <= 0.0 || E->data.Cp0 <= 0.0 ||
       E->data.g0 <= 0.0 || E->data.alpha0 <= 0.0) {
        PyErr_SetString(PyExc_ValueError,
                        "ks, k0, rho0, Cp0, g0, and alpha0 must be positive");
        return NULL;
    }

    E->data.kappa0 = E->data.k0 / (E->data.rho0 * E->data.Cp0);
    E->control.reference_conductivity = E->data.ks / E->data.k0;
    if(E->control.requested_reference_conductivity >= 0.0 &&
       fabs(E->control.requested_reference_conductivity
            - E->control.reference_conductivity)
       > 1.0e-6 * fmax(1.0, fabs(E->control.reference_conductivity))) {
        PyErr_SetString(PyExc_ValueError,
                        "reference_conductivity is derived as ks/k0");
        return NULL;
    }
    E->data.therm_cond = E->data.k0;
    E->data.therm_diff = E->data.kappa0;
    E->data.density = E->data.rho0;
    E->data.Cp = E->data.Cp0;
    E->data.grav_acc = E->data.g0;
    E->data.therm_exp = E->data.alpha0;

    if((E->data.Ttop >= 0.0) != (E->data.Tbottom >= 0.0)) {
        PyErr_SetString(PyExc_ValueError, "Ttop and Tbottom must be specified together");
        return NULL;
    }
    if(E->data.Ttop >= 0.0) {
        if(!isfinite(E->data.Ttop) || !isfinite(E->data.Tbottom) ||
           E->data.Tbottom <= E->data.Ttop) {
            PyErr_SetString(PyExc_ValueError,
                            "Ttop and Tbottom must be finite with Tbottom > Ttop >= 0");
            return NULL;
        }
        if(E->data.ref_temperature > 0.0)
            fprintf(stderr, "WARNING: deltaT is deprecated and ignored; using Tbottom-Ttop\n");
        if(E->control.surface_temp >= 0.0)
            fprintf(stderr, "WARNING: surfaceT is deprecated and ignored; using Ttop/DeltaT\n");
        E->data.ref_temperature = E->data.Tbottom - E->data.Ttop;
        E->control.surface_temp = E->data.Ttop / E->data.ref_temperature;
    }
    else {
        fprintf(stderr, "WARNING: Ttop/Tbottom absent; deltaT/surfaceT compatibility mode is deprecated\n");
        if(E->control.surface_temp < 0.0)
            E->control.surface_temp = 0.1;
    }

    if(E->data.ref_temperature > 0.0)
        E->control.Atemp = E->data.rho0 * E->data.g0 * E->data.alpha0
            * E->data.ref_temperature * radius * radius * radius
            / (E->data.ref_viscosity * E->data.kappa0);
    else
        E->data.ref_temperature = E->control.Atemp * E->data.kappa0
            * E->data.ref_viscosity / (radius * radius * radius)
            / (E->data.rho0 * E->data.g0 * E->data.alpha0);

    if(E->data.Ttop < 0.0) {
        E->data.Ttop = E->control.surface_temp * E->data.ref_temperature;
        E->data.Tbottom = E->data.Ttop + E->data.ref_temperature;
    }
    fprintf(stderr,
            "Temperature reference:\nTtop = %.9g K\nTbottom = %.9g K\nDeltaT = %.9g K\n",
            E->data.Ttop, E->data.Tbottom, E->data.ref_temperature);

    getFloatProperty(properties, "z_lith", E->viscosity.zlith, fp);
    getFloatProperty(properties, "z_410", E->viscosity.z410, fp);
    getFloatProperty(properties, "z_lmantle", E->viscosity.zlm, fp);
    getFloatProperty(properties, "z_cmb", E->viscosity.zcmb, fp); /* this is used as the D" phase change depth */

    /* convert meter to kilometer */
    E->data.radius_km = radius / 1e3;

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;

}



char pyCitcom_IC_set_properties__doc__[] = "";
char pyCitcom_IC_set_properties__name__[] = "IC_set_properties";

PyObject * pyCitcom_IC_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "OOO:IC_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.ic]\n"));

    getIntProperty(properties, "restart", E->control.restart, fp);
    getIntProperty(properties, "post_p", E->control.post_p, fp);
    getIntProperty(properties, "solution_cycles_init", E->monitor.solution_cycles_init, fp);
    getIntProperty(properties, "zero_elapsed_time", E->control.zero_elapsed_time, fp);

    getIntProperty(properties, "tic_method", E->convection.tic_method, fp);

    if (E->convection.tic_method == 0 || E->convection.tic_method == 3) {
	int num_perturb;

	getIntProperty(properties, "num_perturbations", num_perturb, fp);
	if(num_perturb > PERTURB_MAX_LAYERS) {
	    fprintf(stderr, "'num_perturb' greater than allowed value, set to %d\n", PERTURB_MAX_LAYERS);
	    num_perturb = PERTURB_MAX_LAYERS;
	}
	E->convection.number_of_perturbations = num_perturb;

	getIntVectorProperty(properties, "perturbl", E->convection.perturb_ll,
                             num_perturb, fp);
	getIntVectorProperty(properties, "perturbm", E->convection.perturb_mm,
                             num_perturb, fp);
	getIntVectorProperty(properties, "perturblayer", E->convection.load_depth,
                             num_perturb, fp);
	getFloatVectorProperty(properties, "perturbmag", E->convection.perturb_mag,
                               num_perturb, fp);
    }
    else if (E->convection.tic_method == 1) {
        getFloatProperty(properties, "half_space_age", E->convection.half_space_age, fp);
    }
    else if (E->convection.tic_method == 2) {
        getFloatProperty(properties, "half_space_age", E->convection.half_space_age, fp);
        getFloatVectorProperty(properties, "blob_center", E->convection.blob_center, 3, fp);
        if( E->convection.blob_center[0] == -999.0 && E->convection.blob_center[1] == -999.0 && E->convection.blob_center[2] == -999.0 ) {
            E->convection.blob_center[0] = 0.5*(E->control.theta_min+E->control.theta_max);
            E->convection.blob_center[1] = 0.5*(E->control.fi_min+E->control.fi_max);
            E->convection.blob_center[2] = 0.5*(E->sphere.ri+E->sphere.ro);
        }
        getFloatProperty(properties, "blob_radius", E->convection.blob_radius, fp);
        getFloatProperty(properties, "blob_dT", E->convection.blob_dT, fp);
    }

    PUTS(("\n"));

    if (PyErr_Occurred())
      return NULL;

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_Output_set_properties__doc__[] = "";
char pyCitcom_Output_set_properties__name__[] = "Output_set_properties";

PyObject * pyCitcom_Output_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "OOO:Output_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.output]\n"));

    getStringProperty(properties, "output_format", E->output.format, fp);
    getStringProperty(properties, "output_optional", E->output.optional, fp);
    getStringProperty(properties, "profile_optional", E->output.profile_optional, fp);

    getIntProperty(properties, "gzdir_vtkio", E->output.gzdir.vtk_io, fp);
    getIntProperty(properties, "gzdir_rnr", E->output.gzdir.rnr, fp);
    E->output.gzdir.vtk_base_init = 0;
    /* should we save the basis vectors? (memory!) */
    E->output.gzdir.vtk_base_save = 1;

    getIntProperty(properties, "output_ll_max", E->output.llmax, fp);

    getIntProperty(properties, "cb_block_size", E->output.cb_block_size, fp);
    getIntProperty(properties, "cb_buffer_size", E->output.cb_buffer_size, fp);

    getIntProperty(properties, "sieve_buf_size", E->output.sieve_buf_size, fp);

    getIntProperty(properties, "output_alignment", E->output.alignment, fp);
    getIntProperty(properties, "output_alignment_threshold", E->output.alignment_threshold, fp);

    getIntProperty(properties, "cache_mdc_nelmts", E->output.cache_mdc_nelmts, fp);
    getIntProperty(properties, "cache_rdcc_nelmts", E->output.cache_rdcc_nelmts, fp);
    getIntProperty(properties, "cache_rdcc_nbytes", E->output.cache_rdcc_nbytes, fp);

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;

}



char pyCitcom_Param_set_properties__doc__[] = "";
char pyCitcom_Param_set_properties__name__[] = "Param_set_properties";

PyObject * pyCitcom_Param_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "OOO:Param_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.param]\n"));

    getIntProperty(properties, "reference_state", E->refstate.choice, fp);
    if(E->refstate.choice == 0) {
        getStringProperty(properties, "refstate_file", E->refstate.filename, fp);
        getStringProperty(properties, "ala_beta_interval_file",
                          E->refstate.beta_interval_filename, fp);
    }

    getIntProperty(properties, "file_vbcs", E->control.vbcs_file, fp);
    getStringProperty(properties, "vel_bound_file", E->control.velocity_boundary_file, fp);

    getIntProperty(properties, "mat_control", E->control.mat_control, fp);
    getStringProperty(properties, "mat_file", E->control.mat_file, fp);

    getIntProperty(properties, "lith_age", E->control.lith_age, fp);
    getStringProperty(properties, "lith_age_file", E->control.lith_age_file, fp);
    getStringProperty(properties, "flag_depth_new_file", E->control.flag_depth_new_file, fp);
    getStringProperty(properties, "flag_depth_file", E->control.flag_depth_file, fp);
    getStringProperty(properties, "tf_file", E->control.tf_file, fp);
    getIntProperty(properties, "lith_age_time", E->control.lith_age_time, fp);
    getFloatProperty(properties, "lith_age_depth", E->control.lith_age_depth, fp);
    getFloatProperty(properties, "max_plate_age_Ma",
                     E->control.max_plate_age_Ma, fp);
    getFloatProperty(properties, "mantle_temp", E->control.lith_age_mantle_temp, fp);
    getFloatProperty(properties, "bottom_tbl_thickness",
                     E->control.bottom_tbl_thickness, fp);
    getFloatProperty(properties, "bottom_tbl_diffusivity_ratio",
                     E->control.bottom_tbl_diffusivity_ratio, fp);

    getFloatProperty(properties, "start_age", E->control.start_age, fp);
    getIntProperty(properties, "reset_startage", E->control.reset_startage, fp);

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;

}



char pyCitcom_Phase_set_properties__doc__[] = "";
char pyCitcom_Phase_set_properties__name__[] = "Phase_set_properties";

PyObject * pyCitcom_Phase_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    float depth[PHASE_TRANSITIONS];
    float density_jump[PHASE_TRANSITIONS];
    float Ra[PHASE_TRANSITIONS];
    float width[PHASE_TRANSITIONS];
    float clapeyron[PHASE_TRANSITIONS];
    float transT[PHASE_TRANSITIONS];
    int i;

    if (!PyArg_ParseTuple(args, "OOO:Phase_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.phase]\n"));

    getFloatVectorProperty(properties, "phase_depth", depth,
                           PHASE_TRANSITIONS, fp);
    getFloatVectorProperty(properties, "phase_delta_rho", density_jump,
                           PHASE_TRANSITIONS, fp);
    getFloatVectorProperty(properties, "phase_Ra", Ra,
                           PHASE_TRANSITIONS, fp);
    getFloatVectorProperty(properties, "phase_width", width,
                           PHASE_TRANSITIONS, fp);
    getFloatVectorProperty(properties, "phase_clapeyron", clapeyron,
                           PHASE_TRANSITIONS, fp);
    getFloatVectorProperty(properties, "phase_transT", transT,
                           PHASE_TRANSITIONS, fp);

    for(i=0; i<PHASE_TRANSITIONS; i++) {
        E->control.phase[i].depth = depth[i];
        E->control.phase[i].density_jump = density_jump[i];
        E->control.phase[i].Ra = Ra[i];
        E->control.phase[i].clapeyron = clapeyron[i];
        E->control.phase[i].transT = transT[i];
        E->control.phase[i].inv_width = (width[i] == 0.0)? 0.0 : 1.0/width[i];
    }

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;

}



char pyCitcom_Solver_set_properties__doc__[] = "";
char pyCitcom_Solver_set_properties__name__[] = "Solver_set_properties";

PyObject * pyCitcom_Solver_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    float tmp;

    if (!PyArg_ParseTuple(args, "OOO:Solver_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver]\n"));

    getStringProperty(properties, "datadir", E->control.data_dir, fp);
    getStringProperty(properties, "datafile", E->control.data_prefix, fp);
    getStringProperty(properties, "datadir_old", E->control.data_dir_old, fp);
    getStringProperty(properties, "datafile_old", E->control.data_prefix_old, fp);
    getStringProperty(properties, "logfile", E->control.log_template, fp);

    getFloatProperty(properties, "rayleigh", E->control.Atemp, fp);
    getFloatProperty(properties, "dissipation_number", E->control.disptn_number, fp);
    getFloatProperty(properties, "gruneisen", tmp, fp);
     /* special case: if tmp==0, set gruneisen as inf */
     if(tmp != 0)
        E->control.inv_gruneisen = 1/tmp;
    else
        E->control.inv_gruneisen = 0;

    getFloatProperty(properties, "surfaceT", E->control.surface_temp, fp);
    /*getFloatProperty(properties, "adiabaticT0", E->control.adiabaticT0, fp);*/
    getFloatProperty(properties, "Q0", E->control.Q0, fp);

    getIntProperty(properties, "cbf_output_shflux", E->output.cbf_output_shflux, fp);
    getIntProperty(properties, "cbf_output_bhflux", E->output.cbf_output_bhflux, fp);
    getIntProperty(properties, "cbf_use_advection", E->output.cbf_use_advection, fp);

    getIntProperty(properties, "stokes_flow_only", E->control.stokes, fp);

    getIntProperty(properties, "verbose", E->control.verbose, fp);
    getIntProperty(properties, "see_convergence", E->control.print_convergence, fp);

    /* parameters not used in pyre version,
       assigned value here to prevent uninitialized access */
    E->advection.min_timesteps = 1;
    E->advection.max_timesteps = 1;
    E->advection.max_total_timesteps = 1;
    E->control.checkpoint_frequency = 1;
    E->control.record_every = 1;
    E->control.record_all_until = 1;

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_Sphere_set_properties__doc__[] = "";
char pyCitcom_Sphere_set_properties__name__[] = "Sphere_set_properties";

PyObject * pyCitcom_Sphere_set_properties(PyObject *self, PyObject *args)
{
    void full_set_3dsphere_defaults2(struct All_variables *);
    void regional_set_3dsphere_defaults2(struct All_variables *);

    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    int level_factor;

    if (!PyArg_ParseTuple(args, "OOO:Sphere_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.mesher]\n"));

    getIntProperty(properties, "nproc_surf", E->parallel.nprocxy, fp);

    getIntProperty(properties, "nprocx", E->parallel.nprocx, fp);
    getIntProperty(properties, "nprocy", E->parallel.nprocy, fp);
    getIntProperty(properties, "nprocz", E->parallel.nprocz, fp);

    if (E->parallel.nprocxy == 12)
	if (E->parallel.nprocx != E->parallel.nprocy) {
	    char errmsg[] = "!!!! nprocx must equal to nprocy";
	    PyErr_SetString(PyExc_SyntaxError, errmsg);
	    return NULL;
    }

    getIntProperty(properties, "coor", E->control.coor, fp);
    getFloatVectorProperty(properties, "coor_refine", E->control.coor_refine, 4, fp);
    getStringProperty(properties, "coor_file", E->control.coor_file, fp);

    getIntProperty(properties, "nodex", E->mesh.nox, fp);
    getIntProperty(properties, "nodey", E->mesh.noy, fp);
    getIntProperty(properties, "nodez", E->mesh.noz, fp);
    getIntProperty(properties, "levels", E->mesh.levels, fp);

    if (E->mesh.levels < 1 || E->mesh.levels > MAX_LEVELS) {
	char errmsg[] = "!!!! levels must be between 1 and MAX_LEVELS";
	PyErr_SetString(PyExc_SyntaxError, errmsg);
	return NULL;
    }

    level_factor = (int) pow(2.0, E->mesh.levels - 1);
    if (E->parallel.nprocx < 1 || E->parallel.nprocy < 1 ||
	E->parallel.nprocz < 1 ||
	(E->mesh.nox - 1) % (E->parallel.nprocx * level_factor) ||
	(E->mesh.noy - 1) % (E->parallel.nprocy * level_factor) ||
	(E->mesh.noz - 1) % (E->parallel.nprocz * level_factor)) {
	char errmsg[] = "!!!! mesh elements must be divisible by the processor grid times 2^(levels-1)";
	PyErr_SetString(PyExc_SyntaxError, errmsg);
	return NULL;
    }

    E->mesh.mgunitx = (E->mesh.nox - 1) / E->parallel.nprocx /
	level_factor;
    E->mesh.mgunity = (E->mesh.noy - 1) / E->parallel.nprocy /
	level_factor;
    E->mesh.mgunitz = (E->mesh.noz - 1) / E->parallel.nprocz /
	level_factor;

    if (E->mesh.mgunitx < 1 || E->mesh.mgunity < 1 ||
	E->mesh.mgunitz < 1) {
	char errmsg[] = "!!!! coarsest multigrid level must have at least one local element per direction";
	PyErr_SetString(PyExc_SyntaxError, errmsg);
	return NULL;
    }

    if (E->parallel.nprocxy == 12) {
	if (E->mesh.nox != E->mesh.noy) {
	    char errmsg[] = "!!!! nodex must equal to nodey";
	    PyErr_SetString(PyExc_SyntaxError, errmsg);
	    return NULL;
	}
    }

    getDoubleProperty(properties, "radius_outer", E->sphere.ro, fp);
    getDoubleProperty(properties, "radius_inner", E->sphere.ri, fp);


    if (E->parallel.nprocxy == 12) {
        full_set_3dsphere_defaults2(E);
    }
    else {
	getDoubleProperty(properties, "theta_min", E->control.theta_min, fp);
	getDoubleProperty(properties, "theta_max", E->control.theta_max, fp);
	getDoubleProperty(properties, "fi_min", E->control.fi_min, fp);
	getDoubleProperty(properties, "fi_max", E->control.fi_max, fp);

        regional_set_3dsphere_defaults2(E);
    }

    E->mesh.layer[1] = 1;
    E->mesh.layer[2] = 1;
    E->mesh.layer[3] = 1;

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_Tracer_set_properties__doc__[] = "";
char pyCitcom_Tracer_set_properties__name__[] = "Tracer_set_properties";

PyObject * pyCitcom_Tracer_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    double tmp;
    char message[100];

    if (!PyArg_ParseTuple(args, "OOO:Tracer_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.tracer]\n"));

    getIntProperty(properties, "tracer", E->control.tracer, fp);

    getIntProperty(properties, "tracer_enriched", E->control.tracer_enriched, fp);
    if(E->control.tracer_enriched) {
        if(!E->control.tracer)
            myerror(E,"need to switch on tracers for tracer_enriched");

        getFloatProperty(properties, "Q0_enriched", E->control.Q0ER, fp);
        snprintf(message,100,"using compositionally enriched heating: C = 0: %g C = 1: %g (only one composition!)",
                 E->control.Q0,E->control.Q0ER);
        report(E,message);
    }

    getIntProperty(properties, "tracer_ic_method",
                   E->trace.ic_method, fp);

    if (E->trace.ic_method==0) {
        getIntProperty(properties, "tracers_per_element",
                       E->trace.itperel, fp);
    }
    else if (E->trace.ic_method==1) {
        getStringProperty(properties, "tracer_file",
                          E->trace.tracer_file, fp);
    }
    else if (E->trace.ic_method==2) {
    }
    else {
        fprintf(stderr,"Sorry, tracer_ic_method only 0, 1 and 2 available\n");
        fflush(stderr);
        parallel_process_termination();
    }

    getIntProperty(properties, "tracer_flavors", E->trace.nflavors, fp);

    getIntProperty(properties, "ic_method_for_flavors", E->trace.ic_method_for_flavors, fp);
    getIntProperty(properties, "tracer_reclassify_flavors",
                   E->trace.reclassify_flavors, fp);

    if (E->trace.nflavors > 1) {
        switch(E->trace.ic_method_for_flavors){
        case 0:			/* layer */
            E->trace.z_interface = (double*) malloc((E->trace.nflavors-1)
                                                    *sizeof(double));

            getDoubleVectorProperty(properties, "z_interface", E->trace.z_interface, E->trace.nflavors-1, fp);
            break;
        case 1:			/* from grid in top n materials */
            /* file from which to read */
            getStringProperty(properties, "ictracer_grd_file", E->trace.ggrd_file, fp);
            /* which top layers to use */
            getIntProperty(properties, "ictracer_grd_layers", E->trace.ggrd_layers, fp);
            break;
        default:
            fprintf(stderr,"ic_method_for_flavors %i undefined\n",E->trace.ic_method_for_flavors);
            parallel_process_termination();
            break;
        }
    }

    getIntProperty(properties, "itracer_warnings", E->trace.itracer_warnings, fp);

    getIntProperty(properties, "chemical_buoyancy",
                   E->composition.ichemical_buoyancy, fp);

    if (E->control.tracer && E->composition.ichemical_buoyancy==1) {
        getIntProperty(properties, "buoy_type", E->composition.ibuoy_type, fp);

        if (E->composition.ibuoy_type==0)
            E->composition.ncomp = E->trace.nflavors;
        else if (E->composition.ibuoy_type==1)
            E->composition.ncomp = E->trace.nflavors - 1;

        E->composition.buoyancy_ratio = (double*) malloc(E->composition.ncomp
                                                         *sizeof(double));

        getDoubleVectorProperty(properties, "buoyancy_ratio", E->composition.buoyancy_ratio, E->composition.ncomp, fp);
    }


    if(E->parallel.nprocxy == 12) {

        getDoubleProperty(properties, "regular_grid_deltheta", tmp, fp);
        E->trace.deltheta[0] = tmp;
        getDoubleProperty(properties, "regular_grid_delphi", tmp, fp);
        E->trace.delphi[0] = tmp;

        E->trace.ianalytical_tracer_test = 0;
        //getIntProperty(properties, "analytical_tracer_test", E->trace.ianalytical_tracer_test, fp);


        E->composition.icompositional_rheology = 0;
        /*
        getIntProperty(properties, "compositional_rheology", E->composition.icompositional_rheology, fp);

        if (E->composition.icompositional_rheology==1) {
            getDoubleProperty(properties, "compositional_prefactor", E->composition.compositional_rheology_prefactor, fp);
        }
        */
    }
    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}



char pyCitcom_Visc_set_properties__doc__[] = "";
char pyCitcom_Visc_set_properties__name__[] = "Visc_set_properties";

PyObject * pyCitcom_Visc_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;
    int num_mat, i;

    if (!PyArg_ParseTuple(args, "OOO:Visc_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.visc]\n"));

    getStringProperty(properties, "Viscosity", E->viscosity.STRUCTURE, fp);
    if (strcmp(E->viscosity.STRUCTURE,"system") == 0)
	E->viscosity.FROM_SYSTEM = 1;
    else
	E->viscosity.FROM_SYSTEM = 0;

    getIntProperty(properties, "visc_smooth_method", E->viscosity.smooth_cycles, fp);
    getIntProperty(properties, "VISC_UPDATE", E->viscosity.update_allowed, fp);

#define MAX_MAT 40

    getIntProperty(properties, "num_mat", num_mat, fp);
    if(num_mat > MAX_MAT) {
	/* max. allowed material types = 40 */
	fprintf(stderr, "'num_mat' greater than allowed value, set to %d\n", MAX_MAT);
	num_mat = MAX_MAT;
    }
    E->viscosity.num_mat = num_mat;

    getFloatVectorProperty(properties, "visc0",
                           E->viscosity.N0, num_mat, fp);

    getIntProperty(properties, "TDEPV", E->viscosity.TDEPV, fp);
    getIntProperty(properties, "rheol", E->viscosity.RHEOL, fp);
    getFloatVectorProperty(properties, "viscE",
                           E->viscosity.E, num_mat, fp);
    getFloatVectorProperty(properties, "viscT",
                           E->viscosity.T, num_mat, fp);
    getFloatVectorProperty(properties, "viscZ",
                           E->viscosity.Z, num_mat, fp);

    getIntProperty(properties, "SDEPV", E->viscosity.SDEPV, fp);
    getFloatVectorProperty(properties, "sdepv_expt",
                           E->viscosity.sdepv_expt, num_mat, fp);

    getIntProperty(properties, "PDEPV", E->viscosity.PDEPV, fp);
    if (E->viscosity.PDEPV) {
        E->viscosity.pdepv_visited = 0;
        getIntProperty(properties, "pdepv_eff", E->viscosity.pdepv_eff, fp);
        getFloatVectorProperty(properties, "pdepv_a",
                               E->viscosity.pdepv_a, num_mat, fp);
        getFloatVectorProperty(properties, "pdepv_b",
                               E->viscosity.pdepv_b, num_mat, fp);
        getFloatVectorProperty(properties, "pdepv_y",
                               E->viscosity.pdepv_y, num_mat, fp);
        getFloatProperty(properties, "pdepv_offset", E->viscosity.pdepv_offset, fp);
    }
    if(E->viscosity.PDEPV || E->viscosity.SDEPV)
        getFloatProperty(properties, "sdepv_misfit", E->viscosity.sdepv_misfit, fp);

    getIntProperty(properties, "CDEPV", E->viscosity.CDEPV, fp);
    if(E->viscosity.CDEPV){	/* compositional viscosity */
        if(!E->control.tracer)
            myerror(E,"error: CDEPV requires tracers, but tracer is off");
        if(E->trace.nflavors > 30)
            myerror(E,"error: too many flavors for CDEPV");
        /* read in flavor factors */
        getFloatVectorProperty(properties, "cdepv_ff",
                               E->viscosity.cdepv_ff, E->trace.nflavors, fp);
        /* and take the log because we're using a geometric avg */
        for(i=0;i<E->trace.nflavors;i++)
            E->viscosity.cdepv_ff[i] = log(E->viscosity.cdepv_ff[i]);
        E->composition.icompositional_rheology = E->viscosity.CDEPV;
    }

    getIntProperty(properties, "low_visc_channel", E->viscosity.channel, fp);
    getIntProperty(properties, "low_visc_wedge", E->viscosity.wedge, fp);

    getFloatProperty(properties, "lv_min_radius", E->viscosity.lv_min_radius, fp);
    getFloatProperty(properties, "lv_max_radius", E->viscosity.lv_max_radius, fp);
    getFloatProperty(properties, "lv_channel_thickness", E->viscosity.lv_channel_thickness, fp);
    getFloatProperty(properties, "lv_reduction", E->viscosity.lv_reduction, fp);

    getIntProperty(properties, "VMIN", E->viscosity.MIN, fp);
    getFloatProperty(properties, "visc_min", E->viscosity.min_value, fp);

    getIntProperty(properties, "VMAX", E->viscosity.MAX, fp);
    getFloatProperty(properties, "visc_max", E->viscosity.max_value, fp);

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}


char pyCitcom_Incompressible_set_properties__doc__[] = "";
char pyCitcom_Incompressible_set_properties__name__[] = "Incompressible_set_properties";

PyObject * pyCitcom_Incompressible_set_properties(PyObject *self, PyObject *args)
{
    PyObject *obj, *properties, *out;
    struct All_variables *E;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "OOO:Incompressible_set_properties",
			  &obj, &properties, &out))
        return NULL;

    E = (struct All_variables*)(PyCObject_AsVoidPtr(obj));
    fp = get_output_stream(out, E);

    PUTS(("[CitcomS.solver.vsolver]\n"));

    getStringProperty(properties, "Solver", E->control.SOLVER_TYPE, fp);
    getIntProperty(properties, "node_assemble", E->control.NASSEMBLE, fp);
    getIntProperty(properties, "precond", E->control.precondition, fp);

    getDoubleProperty(properties, "accuracy", E->control.accuracy, fp);
    getFloatProperty(properties, "tole_compressibility", E->control.tole_comp, fp);

    getIntProperty(properties, "mg_cycle", E->control.mg_cycle, fp);
    getIntProperty(properties, "down_heavy", E->control.down_heavy, fp);
    getIntProperty(properties, "up_heavy", E->control.up_heavy, fp);

    getIntProperty(properties, "vlowstep", E->control.v_steps_low, fp);
    getIntProperty(properties, "vhighstep", E->control.v_steps_high, fp);
    getIntProperty(properties, "piterations", E->control.p_iterations, fp);

    getIntProperty(properties, "aug_lagr", E->control.augmented_Lagr, fp);
    getDoubleProperty(properties, "aug_number", E->control.augmented, fp);

    getIntProperty(properties, "remove_rigid_rotation", E->control.remove_rigid_rotation, fp);

    getStringProperty(properties, "compressible_formulation",
                      E->control.compressible_formulation, fp);
    if(strcmp(E->control.compressible_formulation, "tala") == 0) {
        E->control.ala_pressure_buoyancy = 0;
    }
    else if(strcmp(E->control.compressible_formulation, "ala") == 0) {
        E->control.ala_pressure_buoyancy = 1;
        if(E->control.inv_gruneisen == 0)
            myerror(E, "compressible_formulation=ala requires gruneisen != 0");
    }
    else {
        myerror(E, "compressible_formulation must be tala or ala");
    }

    getIntProperty(properties, "ala_schur_symmetry_check",
                   E->control.ala_schur_symmetry_check, fp);
    getDoubleProperty(properties, "ala_schur_symmetry_tolerance",
                      E->control.ala_schur_symmetry_tolerance, fp);
    getDoubleProperty(properties, "ala_augmented_lagrangian_gamma",
                      E->control.ala_augmented_lagrangian_gamma, fp);
    getStringProperty(properties, "ala_beta_element_source",
                      E->control.ala_beta_element_source, fp);
    getDoubleProperty(properties, "ala_inner_accuracy_max",
                      E->control.ala_inner_accuracy_max, fp);
    getDoubleProperty(properties, "ala_inner_accuracy_factor",
                      E->control.ala_inner_accuracy_factor, fp);
    getIntProperty(properties, "ala_pcg_restart_interval",
                   E->control.ala_pcg_restart_interval, fp);
    getStringProperty(properties, "ala_outer_solver",
                      E->control.ala_outer_solver, fp);
    getDoubleProperty(
        properties, "ala_coupled_initial_velocity_relative_tolerance",
        E->control.ala_coupled_initial_velocity_relative_tolerance, fp);
    getDoubleProperty(properties, "ala_coupled_inner_relative_tolerance",
                      E->control.ala_coupled_inner_relative_tolerance, fp);
    getIntProperty(properties, "ala_coupled_inner_max_cycles",
                   E->control.ala_coupled_inner_max_cycles, fp);
    getIntProperty(properties, "ala_coupled_inner_progress_interval",
                   E->control.ala_coupled_inner_progress_interval, fp);
    getIntProperty(properties, "ala_coupled_defect_corrections",
                   E->control.ala_coupled_defect_corrections, fp);
    getIntProperty(properties, "ala_coupled_multilevel_audit_only",
                   E->control.ala_coupled_multilevel_audit_only, fp);
    getIntProperty(properties, "ala_coupled_first_preconditioner_audit_only",
                   E->control.ala_coupled_first_preconditioner_audit_only, fp);
    getIntProperty(properties, "ala_coupled_debug_stop_iteration",
                   E->control.ala_coupled_debug_stop_iteration, fp);
    getIntProperty(properties, "ala_coupled_element_vanka",
                   E->control.ala_coupled_element_vanka, fp);
    getIntProperty(properties, "ala_coupled_multilevel_vcycle",
                   E->control.ala_coupled_multilevel_vcycle, fp);
    getIntProperty(properties, "ala_coupled_multilevel_coarse_sweeps",
                   E->control.ala_coupled_multilevel_coarse_sweeps, fp);
    getDoubleProperty(properties, "ala_coupled_multilevel_coarse_weight",
                      E->control.ala_coupled_multilevel_coarse_weight, fp);
    getIntProperty(properties, "ala_viscosity_spectrum_diagnostics",
                   E->control.ala_viscosity_spectrum_diagnostics, fp);
    getIntProperty(properties, "ala_viscosity_spectrum_interval",
                   E->control.ala_viscosity_spectrum_interval, fp);
    getIntProperty(properties, "ala_coupled_shallow_vanka_layers",
                   E->control.ala_coupled_shallow_vanka_layers, fp);
    getIntProperty(properties, "ala_coupled_shallow_vanka_core_layers",
                   E->control.ala_coupled_shallow_vanka_core_layers, fp);
    getIntProperty(properties, "ala_coupled_shallow_vanka_band_sweeps",
                   E->control.ala_coupled_shallow_vanka_band_sweeps, fp);
    getIntProperty(properties, "ala_coupled_shallow_vanka_sweeps",
                   E->control.ala_coupled_shallow_vanka_sweeps, fp);
    getIntProperty(properties, "ala_coupled_shallow_vanka_warm_sweeps",
                   E->control.ala_coupled_shallow_vanka_warm_sweeps, fp);
    getDoubleProperty(properties, "ala_unaugmented_momentum_tolerance",
                      E->control.ala_unaugmented_momentum_tolerance, fp);
    getIntProperty(properties, "ala_feasibility_audit",
                   E->control.ala_feasibility_audit, fp);
    getIntProperty(properties, "ala_feasibility_window",
                   E->control.ala_feasibility_window, fp);
    getDoubleProperty(properties, "ala_feasibility_min_reduction",
                      E->control.ala_feasibility_min_reduction, fp);
    getIntProperty(properties, "ala_hybrid_convergence",
                   E->control.ala_hybrid_convergence, fp);
    getDoubleProperty(properties, "ala_div_v_tolerance",
                      E->control.ala_div_v_tolerance, fp);
    getDoubleProperty(properties, "ala_update_tolerance",
                      E->control.ala_update_tolerance, fp);
    getIntProperty(properties, "ala_consecutive_steps",
                   E->control.ala_consecutive_steps, fp);
    getIntProperty(properties, "ala_depth_diagnostics",
                   E->control.ala_depth_diagnostics, fp);
    getIntProperty(properties, "ala_depth_diagnostic_interval",
                   E->control.ala_depth_diagnostic_interval, fp);
    getIntProperty(properties, "ala_depth_diagnostic_bins",
                   E->control.ala_depth_diagnostic_bins, fp);
    getIntProperty(properties, "ala_beta_causal_diagnostics",
                   E->control.ala_beta_causal_diagnostics, fp);
    getIntProperty(properties, "ala_coarse_residual_diagnostics",
                   E->control.ala_coarse_residual_diagnostics, fp);
    getIntProperty(properties, "ala_coarse_residual_interval",
                   E->control.ala_coarse_residual_interval, fp);
    getIntProperty(properties, "ala_coarse_residual_levels",
                   E->control.ala_coarse_residual_levels, fp);
    getIntProperty(properties, "ala_two_level_preconditioner",
                   E->control.ala_two_level_preconditioner, fp);
    getIntProperty(properties, "ala_two_level_offset",
                   E->control.ala_two_level_offset, fp);
    getIntProperty(properties, "ala_two_level_coarse_iterations",
                   E->control.ala_two_level_coarse_iterations, fp);
    getDoubleProperty(properties, "ala_two_level_coarse_damping",
                      E->control.ala_two_level_coarse_damping, fp);
    getStringProperty(properties, "ala_two_level_coarse_solver",
                      E->control.ala_two_level_coarse_solver, fp);
    getDoubleProperty(properties, "ala_two_level_coarse_eigenvalue_min",
                      E->control.ala_two_level_coarse_eigenvalue_min, fp);
    getDoubleProperty(properties, "ala_two_level_coarse_eigenvalue_max",
                      E->control.ala_two_level_coarse_eigenvalue_max, fp);
    getDoubleProperty(properties, "ala_two_level_coarse_weight",
                      E->control.ala_two_level_coarse_weight, fp);
    getStringProperty(properties, "ala_two_level_velocity_solver",
                      E->control.ala_two_level_velocity_solver, fp);
    getIntProperty(properties, "ala_two_level_velocity_iterations",
                   E->control.ala_two_level_velocity_iterations, fp);
    getDoubleProperty(properties, "ala_two_level_velocity_eigenvalue_min",
                      E->control.ala_two_level_velocity_eigenvalue_min, fp);
    getDoubleProperty(properties, "ala_two_level_velocity_eigenvalue_max",
                      E->control.ala_two_level_velocity_eigenvalue_max, fp);
    getIntProperty(properties, "ala_pressure_multigrid",
                   E->control.ala_pressure_multigrid, fp);
    getIntProperty(properties, "ala_pressure_multigrid_min_level",
                   E->control.ala_pressure_multigrid_min_level, fp);
    getIntProperty(properties, "ala_pressure_multigrid_pre_smooth",
                   E->control.ala_pressure_multigrid_pre_smooth, fp);
    getIntProperty(properties, "ala_pressure_multigrid_post_smooth",
                   E->control.ala_pressure_multigrid_post_smooth, fp);
    getIntProperty(properties, "ala_pressure_multigrid_coarse_iterations",
                   E->control.ala_pressure_multigrid_coarse_iterations, fp);
    getDoubleProperty(properties, "ala_pressure_multigrid_damping",
                      E->control.ala_pressure_multigrid_damping, fp);
    getDoubleProperty(properties, "ala_pressure_multigrid_weight",
                      E->control.ala_pressure_multigrid_weight, fp);
    getIntProperty(properties, "ala_global_coarse_preconditioner",
                   E->control.ala_global_coarse_preconditioner, fp);
    getDoubleProperty(properties, "ala_global_coarse_weight",
                      E->control.ala_global_coarse_weight, fp);
    getDoubleProperty(properties, "ala_global_coarse_regularization",
                      E->control.ala_global_coarse_regularization, fp);
    getIntProperty(properties, "ala_shallow_patch_preconditioner",
                   E->control.ala_shallow_patch_preconditioner, fp);
    getDoubleProperty(properties, "ala_shallow_patch_depth_km",
                      E->control.ala_shallow_patch_depth_km, fp);
    getDoubleProperty(properties, "ala_shallow_patch_weight",
                      E->control.ala_shallow_patch_weight, fp);
    getDoubleProperty(properties, "ala_shallow_patch_regularization",
                      E->control.ala_shallow_patch_regularization, fp);
    getIntProperty(properties, "ala_shallow_patch_horizontal_elements",
                   E->control.ala_shallow_patch_horizontal_elements, fp);
    getIntProperty(properties, "ala_shallow_patch_horizontal_stride",
                   E->control.ala_shallow_patch_horizontal_stride, fp);
    getIntProperty(properties, "ala_shallow_patch_mpi_overlap",
                   E->control.ala_shallow_patch_mpi_overlap, fp);
    getStringProperty(properties, "ala_shallow_patch_velocity_solver",
                      E->control.ala_shallow_patch_velocity_solver, fp);
    getIntProperty(properties, "ala_geneo_preconditioner",
                   E->control.ala_geneo_preconditioner, fp);
    getStringProperty(properties, "ala_geneo_basis_type",
                      E->control.ala_geneo_basis_type, fp);
    getDoubleProperty(properties, "ala_geneo_eigenvalue_threshold",
                      E->control.ala_geneo_eigenvalue_threshold, fp);
    getIntProperty(properties, "ala_geneo_min_modes_per_rank",
                   E->control.ala_geneo_min_modes_per_rank, fp);
    getIntProperty(properties, "ala_geneo_max_modes_per_rank",
                   E->control.ala_geneo_max_modes_per_rank, fp);
    getIntProperty(properties, "ala_geneo_horizontal_bins",
                   E->control.ala_geneo_horizontal_bins, fp);
    getIntProperty(properties, "ala_geneo_radial_bins",
                   E->control.ala_geneo_radial_bins, fp);
    getIntProperty(properties, "ala_geneo_rank_group_x",
                   E->control.ala_geneo_rank_group_x, fp);
    getIntProperty(properties, "ala_geneo_rank_group_y",
                   E->control.ala_geneo_rank_group_y, fp);
    getIntProperty(properties, "ala_geneo_max_global_modes",
                   E->control.ala_geneo_max_global_modes, fp);
    getDoubleProperty(properties, "ala_geneo_weight",
                      E->control.ala_geneo_weight, fp);
    getDoubleProperty(properties, "ala_geneo_regularization",
                      E->control.ala_geneo_regularization, fp);
    getIntProperty(properties, "ala_radial_line_preconditioner",
                   E->control.ala_radial_line_preconditioner, fp);
    getIntProperty(properties, "ala_element_vanka_smoother",
                   E->control.ala_element_vanka_smoother, fp);
    getDoubleProperty(properties, "ala_element_vanka_damping",
                      E->control.ala_element_vanka_damping, fp);
    getDoubleProperty(properties, "ala_element_vanka_pressure_damping",
                      E->control.ala_element_vanka_pressure_damping, fp);
    getDoubleProperty(properties, "ala_element_vanka_regularization",
                      E->control.ala_element_vanka_regularization, fp);
    getIntProperty(properties, "ala_element_vanka_galerkin_schur",
                   E->control.ala_element_vanka_galerkin_schur, fp);
    getIntProperty(properties, "ala_element_vanka_rebuild_interval",
                   E->control.ala_element_vanka_rebuild_interval, fp);
    if(E->control.ala_schur_symmetry_tolerance <= 0.0)
        myerror(E, "ala_schur_symmetry_tolerance must be positive");
    if(E->control.ala_augmented_lagrangian_gamma < 0.0)
        myerror(E, "ala_augmented_lagrangian_gamma must be nonnegative");
    if(E->control.ala_augmented_lagrangian_gamma > 0.0 &&
       !E->control.ala_pressure_buoyancy)
        myerror(E, "ala_augmented_lagrangian_gamma requires "
                "compressible_formulation=ala");
    if(E->control.ala_augmented_lagrangian_gamma > 0.0 &&
       E->control.augmented_Lagr)
        myerror(E, "ala_augmented_lagrangian_gamma and aug_lagr are "
                "mutually exclusive");
    if(strcmp(E->control.ala_beta_element_source,"supplied_average") != 0 &&
       strcmp(E->control.ala_beta_element_source,"density_log_secant") != 0 &&
       strcmp(E->control.ala_beta_element_source,"interval") != 0)
        myerror(E, "ala_beta_element_source must be supplied_average, "
                "density_log_secant, or interval");
    if(strcmp(E->control.ala_beta_element_source,"interval") == 0 &&
       (!E->control.ala_pressure_buoyancy || E->refstate.choice != 0))
        myerror(E, "ala_beta_element_source=interval requires "
                "compressible_formulation=ala and reference_state=0");
    if(E->control.ala_inner_accuracy_max <= 0.0)
        myerror(E, "ala_inner_accuracy_max must be positive");
    if(E->control.ala_inner_accuracy_factor <= 0.0)
        myerror(E, "ala_inner_accuracy_factor must be positive");
    if(E->control.ala_pcg_restart_interval < 1)
        myerror(E, "ala_pcg_restart_interval must be at least one");
    if(E->control.ala_coupled_initial_velocity_relative_tolerance < 0.0 ||
       E->control.ala_coupled_initial_velocity_relative_tolerance >= 1.0)
        myerror(E,
                "ala_coupled_initial_velocity_relative_tolerance must be "
                "in [0,1)");
    if(E->control.ala_coupled_inner_relative_tolerance <= 0.0 ||
       E->control.ala_coupled_inner_relative_tolerance >= 1.0)
        myerror(E, "ala_coupled_inner_relative_tolerance must be in (0,1)");
    if(E->control.ala_coupled_inner_max_cycles < 1)
        myerror(E, "ala_coupled_inner_max_cycles must be at least one");
    if(E->control.ala_coupled_inner_progress_interval < 1)
        myerror(E, "ala_coupled_inner_progress_interval must be at least one");
    if(E->control.ala_coupled_defect_corrections < 0 ||
       E->control.ala_coupled_defect_corrections > 2)
        myerror(E,
                "ala_coupled_defect_corrections must be between zero and two");
    if(strcmp(E->control.ala_outer_solver,"pcg") != 0 &&
       strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_outer_solver must be pcg, fgmres, or coupled_fgmres");
    if(E->control.ala_coupled_defect_corrections > 0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_coupled_defect_corrections requires "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_initial_velocity_relative_tolerance > 0.0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E,
                "ala_coupled_initial_velocity_relative_tolerance requires "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_multilevel_audit_only &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_coupled_multilevel_audit_only requires "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_first_preconditioner_audit_only &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_coupled_first_preconditioner_audit_only requires "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_debug_stop_iteration < 0)
        myerror(E, "ala_coupled_debug_stop_iteration must be nonnegative");
    if(E->control.ala_coupled_debug_stop_iteration > 0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_coupled_debug_stop_iteration requires "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_element_vanka &&
       (!E->control.ala_element_vanka_smoother ||
        strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0))
        myerror(E, "ala_coupled_element_vanka requires "
                "ala_element_vanka_smoother=on and "
                "ala_outer_solver=coupled_fgmres");
    if(E->control.ala_coupled_element_vanka &&
       E->control.ala_coupled_defect_corrections>0)
        myerror(E, "ala_coupled_element_vanka currently requires "
                "ala_coupled_defect_corrections=0");
    if(E->control.ala_coupled_multilevel_vcycle &&
       !E->control.ala_coupled_element_vanka)
        myerror(E, "ala_coupled_multilevel_vcycle requires "
                "ala_coupled_element_vanka=on");
    if(E->control.ala_coupled_multilevel_coarse_sweeps < 2 ||
       E->control.ala_coupled_multilevel_coarse_sweeps > 100)
        myerror(E,
                "ala_coupled_multilevel_coarse_sweeps must be in [2,100]");
    if(E->control.ala_coupled_multilevel_coarse_weight <= 0.0 ||
       E->control.ala_coupled_multilevel_coarse_weight > 1.0)
        myerror(E, "ala_coupled_multilevel_coarse_weight must be in (0,1]");
    if(E->control.ala_viscosity_spectrum_interval < 1)
        myerror(E, "ala_viscosity_spectrum_interval must be at least one");
    if(E->control.ala_coupled_shallow_vanka_layers < 0 ||
       E->control.ala_coupled_shallow_vanka_layers > 64)
        myerror(E, "ala_coupled_shallow_vanka_layers must be in [0,64]");
    if(E->control.ala_coupled_shallow_vanka_core_layers < -1 ||
       E->control.ala_coupled_shallow_vanka_core_layers >
          E->control.ala_coupled_shallow_vanka_layers)
        myerror(E, "ala_coupled_shallow_vanka_core_layers must be -1 or in "
                "[0,ala_coupled_shallow_vanka_layers]");
    if(E->control.ala_coupled_shallow_vanka_band_sweeps < 0 ||
       E->control.ala_coupled_shallow_vanka_band_sweeps > 8)
        myerror(E, "ala_coupled_shallow_vanka_band_sweeps must be in [0,8]");
    if(E->control.ala_coupled_shallow_vanka_band_sweeps > 0 &&
       (E->control.ala_coupled_shallow_vanka_core_layers < 0 ||
        E->control.ala_coupled_shallow_vanka_core_layers >=
          E->control.ala_coupled_shallow_vanka_layers))
        myerror(E, "ala_coupled_shallow_vanka_band_sweeps requires "
                "0 <= core_layers < shallow_layers");
    if(E->control.ala_coupled_shallow_vanka_sweeps < 0 ||
       E->control.ala_coupled_shallow_vanka_sweeps > 32)
        myerror(E, "ala_coupled_shallow_vanka_sweeps must be in [0,32]");
    if(E->control.ala_coupled_shallow_vanka_warm_sweeps < -1 ||
       E->control.ala_coupled_shallow_vanka_warm_sweeps >
          E->control.ala_coupled_shallow_vanka_sweeps)
        myerror(E, "ala_coupled_shallow_vanka_warm_sweeps must be -1 or in "
                "[0,ala_coupled_shallow_vanka_sweeps]");
    if((E->control.ala_coupled_shallow_vanka_layers == 0) !=
       (E->control.ala_coupled_shallow_vanka_sweeps == 0))
        myerror(E, "ala_coupled_shallow_vanka_layers and sweeps must both "
                "be zero or both be positive");
    if(E->control.ala_coupled_shallow_vanka_sweeps > 0 &&
       !E->control.ala_coupled_multilevel_vcycle)
        myerror(E, "ala_coupled_shallow_vanka_sweeps requires "
                "ala_coupled_multilevel_vcycle=on");
    if(E->control.ala_coupled_shallow_vanka_warm_sweeps >= 0 &&
       E->control.ala_coupled_shallow_vanka_layers == 0)
        myerror(E, "ala_coupled_shallow_vanka_warm_sweeps requires "
                "shallow Vanka layers and cold sweeps");
    if(E->control.ala_unaugmented_momentum_tolerance < 0.0)
        myerror(E, "ala_unaugmented_momentum_tolerance must be nonnegative");
    if(E->control.ala_unaugmented_momentum_tolerance > 0.0 &&
       strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
       strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)
        myerror(E, "ala_unaugmented_momentum_tolerance requires "
                "ala_outer_solver=fgmres or coupled_fgmres");
    if(E->control.ala_feasibility_window < 1)
        myerror(E, "ala_feasibility_window must be at least one");
    if(E->control.ala_feasibility_min_reduction < 0.0 ||
       E->control.ala_feasibility_min_reduction >= 1.0)
        myerror(E, "ala_feasibility_min_reduction must be in [0,1)");
    if(E->control.ala_div_v_tolerance <= 0.0)
        myerror(E, "ala_div_v_tolerance must be positive");
    if(E->control.ala_update_tolerance <= 0.0)
        myerror(E, "ala_update_tolerance must be positive");
    if(E->control.ala_consecutive_steps < 1)
        myerror(E, "ala_consecutive_steps must be at least one");
    if(E->control.ala_depth_diagnostic_interval < 1)
        myerror(E, "ala_depth_diagnostic_interval must be at least one");
    if(E->control.ala_depth_diagnostic_bins < 1 ||
       E->control.ala_depth_diagnostic_bins > 128)
        myerror(E, "ala_depth_diagnostic_bins must be between 1 and 128");
    if(E->control.ala_coarse_residual_interval < 1)
        myerror(E, "ala_coarse_residual_interval must be at least one");
    if(E->control.ala_coarse_residual_levels < 1 ||
       E->control.ala_coarse_residual_levels > 10)
        myerror(E, "ala_coarse_residual_levels must be between 1 and 10");
    if(E->control.ala_two_level_offset < 1 ||
       E->control.ala_two_level_offset > 10)
        myerror(E, "ala_two_level_offset must be between 1 and 10");
    if(E->control.ala_two_level_coarse_iterations < 1)
        myerror(E, "ala_two_level_coarse_iterations must be at least one");
    if(E->control.ala_two_level_coarse_damping <= 0.0 ||
       E->control.ala_two_level_coarse_damping >= 2.0/27.0)
        myerror(E, "ala_two_level_coarse_damping must be in (0,2/27)");
    if(strcmp(E->control.ala_two_level_coarse_solver,"jacobi") != 0 &&
       strcmp(E->control.ala_two_level_coarse_solver,"chebyshev") != 0)
        myerror(E, "ala_two_level_coarse_solver must be jacobi or chebyshev");
    if(E->control.ala_two_level_coarse_eigenvalue_min <= 0.0 ||
       E->control.ala_two_level_coarse_eigenvalue_max <=
         E->control.ala_two_level_coarse_eigenvalue_min)
        myerror(E, "ALA two-level Chebyshev eigenvalue interval is invalid");
    if(E->control.ala_two_level_coarse_weight <= 0.0 ||
       E->control.ala_two_level_coarse_weight > 1.0)
        myerror(E, "ala_two_level_coarse_weight must be in (0,1]");
    if(strcmp(E->control.ala_two_level_velocity_solver,"diagonal") != 0 &&
       strcmp(E->control.ala_two_level_velocity_solver,"chebyshev") != 0)
        myerror(E, "ala_two_level_velocity_solver must be diagonal or chebyshev");
    if(E->control.ala_two_level_velocity_iterations < 1)
        myerror(E, "ala_two_level_velocity_iterations must be at least one");
    if(E->control.ala_two_level_velocity_eigenvalue_min <= 0.0 ||
       E->control.ala_two_level_velocity_eigenvalue_max <=
         E->control.ala_two_level_velocity_eigenvalue_min)
        myerror(E, "ALA two-level velocity Chebyshev eigenvalue interval is invalid");
    if(E->control.ala_pressure_multigrid &&
       (E->control.ala_pressure_multigrid_min_level < E->mesh.gridmin ||
        E->control.ala_pressure_multigrid_min_level > E->mesh.levmax))
        myerror(E, "ala_pressure_multigrid_min_level is outside the mesh hierarchy");
    if(E->control.ala_pressure_multigrid &&
       (E->control.ala_pressure_multigrid_pre_smooth < 1 ||
        E->control.ala_pressure_multigrid_post_smooth < 1))
        myerror(E, "ALA pressure multigrid smoothing counts must be positive");
    if(E->control.ala_pressure_multigrid &&
       E->control.ala_pressure_multigrid_coarse_iterations < 1)
        myerror(E, "ALA pressure multigrid coarse iterations must be positive");
    if(E->control.ala_pressure_multigrid &&
       (E->control.ala_pressure_multigrid_damping <= 0.0 ||
        E->control.ala_pressure_multigrid_damping >= 2.0/27.0))
        myerror(E, "ala_pressure_multigrid_damping must be in (0,2/27)");
    if(E->control.ala_pressure_multigrid &&
       (E->control.ala_pressure_multigrid_weight <= 0.0 ||
        E->control.ala_pressure_multigrid_weight > 1.0))
        myerror(E, "ala_pressure_multigrid_weight must be in (0,1]");
    if(E->control.ala_global_coarse_weight <= 0.0 ||
       E->control.ala_global_coarse_weight > 1.0)
        myerror(E, "ala_global_coarse_weight must be in (0,1]");
    if(E->control.ala_global_coarse_regularization < 0.0 ||
       E->control.ala_global_coarse_regularization > 1.0e-4)
        myerror(E, "ala_global_coarse_regularization must be in [0,1e-4]");
    if(E->control.ala_shallow_patch_depth_km <= 0.0)
        myerror(E, "ala_shallow_patch_depth_km must be positive");
    if(E->control.ala_shallow_patch_weight <= 0.0 ||
       E->control.ala_shallow_patch_weight > 1.0)
        myerror(E, "ala_shallow_patch_weight must be in (0,1]");
    if(E->control.ala_shallow_patch_regularization < 0.0 ||
       E->control.ala_shallow_patch_regularization > 0.1)
        myerror(E, "ala_shallow_patch_regularization must be in [0,0.1]");
    if(E->control.ala_shallow_patch_horizontal_elements < 2 ||
       E->control.ala_shallow_patch_horizontal_elements > 8)
        myerror(E, "ala_shallow_patch_horizontal_elements must be in [2,8]");
    if(E->control.ala_shallow_patch_horizontal_stride < 1 ||
       E->control.ala_shallow_patch_horizontal_stride >
         E->control.ala_shallow_patch_horizontal_elements)
        myerror(E, "ala_shallow_patch_horizontal_stride must be in "
                "[1,ala_shallow_patch_horizontal_elements]");
    if(E->control.ala_shallow_patch_mpi_overlap < 1 ||
       E->control.ala_shallow_patch_mpi_overlap > 4)
        myerror(E, "ala_shallow_patch_mpi_overlap must be in [1,4]");
    if(2*E->control.ala_shallow_patch_mpi_overlap >
       E->control.ala_shallow_patch_horizontal_elements)
        myerror(E, "twice ala_shallow_patch_mpi_overlap must not exceed "
                "ala_shallow_patch_horizontal_elements");
    if(strcmp(E->control.ala_shallow_patch_velocity_solver,"diagonal") != 0 &&
       strcmp(E->control.ala_shallow_patch_velocity_solver,"node_block") != 0)
        myerror(E, "ala_shallow_patch_velocity_solver must be diagonal or "
                "node_block");
    if(E->control.ala_geneo_eigenvalue_threshold <= 0.0)
        myerror(E, "ala_geneo_eigenvalue_threshold must be positive");
    if(strcmp(E->control.ala_geneo_basis_type,"spectral") != 0 &&
       strcmp(E->control.ala_geneo_basis_type,"radial_partition") != 0)
        myerror(E, "ala_geneo_basis_type must be spectral or "
                "radial_partition");
    if(strcmp(E->control.ala_geneo_basis_type,"radial_partition") == 0 &&
       (E->control.ala_geneo_min_modes_per_rank
          != E->control.ala_geneo_radial_bins ||
        E->control.ala_geneo_max_modes_per_rank
          != E->control.ala_geneo_radial_bins))
        myerror(E, "radial_partition requires GenEO min/max modes equal "
                "ala_geneo_radial_bins");
    if(strcmp(E->control.ala_geneo_basis_type,"radial_partition") == 0 &&
       E->control.ala_geneo_rank_group_x == 1 &&
       E->control.ala_geneo_rank_group_y == 1)
        myerror(E, "radial_partition requires a cross-rank GenEO group");
    if(E->control.ala_geneo_min_modes_per_rank < 1 ||
       E->control.ala_geneo_max_modes_per_rank <
         E->control.ala_geneo_min_modes_per_rank ||
       E->control.ala_geneo_max_modes_per_rank > 8)
        myerror(E, "ALA GenEO modes per rank must satisfy 1 <= min <= max <= 8");
    if(E->control.ala_geneo_horizontal_bins < 2 ||
       E->control.ala_geneo_horizontal_bins > 8 ||
       E->control.ala_geneo_radial_bins < 1 ||
       E->control.ala_geneo_radial_bins > 4)
        myerror(E, "ALA GenEO bin counts are outside supported bounds");
    if(E->control.ala_geneo_rank_group_x < 1 ||
       E->control.ala_geneo_rank_group_x > E->parallel.nprocx ||
       E->control.ala_geneo_rank_group_y < 1 ||
       E->control.ala_geneo_rank_group_y > E->parallel.nprocy)
        myerror(E,
                "ALA GenEO rank groups must fit the horizontal processor grid");
    if(E->control.ala_geneo_horizontal_bins *
         E->control.ala_geneo_rank_group_x *
         E->control.ala_geneo_horizontal_bins *
         E->control.ala_geneo_rank_group_y *
         E->control.ala_geneo_radial_bins > 256)
        myerror(E,
                "ALA GenEO cross-rank aggregate exceeds 256 spectral bins");
    if(E->control.ala_geneo_max_global_modes < 1 ||
       E->control.ala_geneo_max_global_modes > 4096)
        myerror(E, "ala_geneo_max_global_modes must be between 1 and 4096");
    if(E->control.ala_geneo_weight <= 0.0 ||
       E->control.ala_geneo_weight > 1.0)
        myerror(E, "ala_geneo_weight must be in (0,1]");
    if(E->control.ala_geneo_regularization < 0.0 ||
       E->control.ala_geneo_regularization > 1.0e-3)
        myerror(E, "ala_geneo_regularization must be in [0,1e-3]");
    if(E->control.ala_element_vanka_damping <= 0.0 ||
       E->control.ala_element_vanka_damping > 1.0)
        myerror(E, "ala_element_vanka_damping must be in (0,1]");
    if(E->control.ala_element_vanka_pressure_damping <= 0.0 ||
       E->control.ala_element_vanka_pressure_damping > 1.0)
        myerror(E, "ala_element_vanka_pressure_damping must be in (0,1]");
    if(E->control.ala_element_vanka_regularization < 0.0 ||
       E->control.ala_element_vanka_regularization > 0.1)
        myerror(E, "ala_element_vanka_regularization must be in [0,0.1]");
    if(E->control.ala_element_vanka_rebuild_interval < 1 ||
       E->control.ala_element_vanka_rebuild_interval > 100)
        myerror(E, "ala_element_vanka_rebuild_interval must be in [1,100]");

    if(E->control.inv_gruneisen != 0) {
        /* "cg" is legacy split; strict ALA uses bicg or ala_cg. */
        getStringProperty(properties, "uzawa", E->control.uzawa, fp);
        if(strcmp(E->control.uzawa, "cg") == 0) {
            if(E->control.ala_pressure_buoyancy)
                myerror(E,
                        "strict ALA requires uzawa=bicg or uzawa=ala_cg");
            /* more convergence parameters for "cg" */
            getIntProperty(properties, "compress_iter_maxstep", E->control.compress_iter_maxstep, fp);
            getFloatProperty(properties, "relative_err_accuracy", E->control.relative_err_accuracy, fp);
        }
        else if(strcmp(E->control.uzawa, "bicg") == 0) {
            getIntProperty(properties, "compress_iter_maxstep", E->control.compress_iter_maxstep, fp);
            getFloatProperty(properties, "relative_err_accuracy", E->control.relative_err_accuracy, fp);
        }
        else if(strcmp(E->control.uzawa, "ala_cg") == 0) {
            if(!E->control.ala_pressure_buoyancy)
                myerror(E,
                        "uzawa=ala_cg requires compressible_formulation=ala");
            getIntProperty(properties, "compress_iter_maxstep", E->control.compress_iter_maxstep, fp);
            getFloatProperty(properties, "relative_err_accuracy", E->control.relative_err_accuracy, fp);
        }
        else
            myerror(E, "Error: unknown Uzawa iteration");
    }
    if((strcmp(E->control.ala_outer_solver,"fgmres") == 0 ||
        strcmp(E->control.ala_outer_solver,"coupled_fgmres") == 0) &&
       strcmp(E->control.uzawa,"ala_cg") != 0)
        myerror(E, "ALA FGMRES outer solvers require uzawa=ala_cg");
    if(E->control.ala_radial_line_preconditioner &&
       (!E->control.precondition ||
        !E->control.ala_pressure_buoyancy ||
        strcmp(E->control.uzawa,"ala_cg") != 0))
        myerror(E,
                "ala_radial_line_preconditioner requires precond=on, "
                "compressible_formulation=ala, and uzawa=ala_cg");
    if(E->control.ala_two_level_preconditioner &&
       (!E->control.precondition ||
        !E->control.ala_pressure_buoyancy ||
        strcmp(E->control.uzawa,"ala_cg") != 0))
        myerror(E,
                "ala_two_level_preconditioner requires precond=on, "
                "compressible_formulation=ala, and uzawa=ala_cg");
    if(E->control.ala_two_level_preconditioner &&
       E->control.ala_radial_line_preconditioner)
        myerror(E, "ALA two-level and radial-line preconditioners are "
                "mutually exclusive");
    if(E->control.ala_pressure_multigrid &&
       (!E->control.precondition ||
        !E->control.ala_pressure_buoyancy ||
        E->control.ala_augmented_lagrangian_gamma <= 0.0 ||
        (strcmp(E->control.ala_outer_solver,"fgmres") != 0 &&
         strcmp(E->control.ala_outer_solver,"coupled_fgmres") != 0)))
        myerror(E, "ala_pressure_multigrid requires precond=on, strict ALA, "
                "positive gamma, and an FGMRES outer solver");
    if(E->control.ala_pressure_multigrid &&
       E->control.ala_two_level_preconditioner)
        myerror(E, "ALA pressure multigrid and legacy two-level "
                "preconditioners are mutually exclusive");
    if(E->control.ala_global_coarse_preconditioner &&
       !E->control.ala_two_level_preconditioner)
        myerror(E, "ala_global_coarse_preconditioner requires "
                "ala_two_level_preconditioner=on");
    if(E->control.ala_shallow_patch_preconditioner &&
       (!E->control.precondition ||
        !E->control.ala_pressure_buoyancy ||
        strcmp(E->control.uzawa,"ala_cg") != 0))
        myerror(E,
                "ala_shallow_patch_preconditioner requires precond=on, "
                "compressible_formulation=ala, and uzawa=ala_cg");
    if(E->control.ala_shallow_patch_preconditioner &&
       E->control.ala_radial_line_preconditioner)
        myerror(E, "ALA shallow-patch and radial-line preconditioners are "
                "mutually exclusive");
    if(E->control.ala_geneo_preconditioner &&
       (!E->control.precondition ||
        !E->control.ala_pressure_buoyancy ||
        strcmp(E->control.uzawa,"ala_cg") != 0 ||
        !E->control.ala_shallow_patch_preconditioner))
        myerror(E, "ala_geneo_preconditioner requires strict ALA PCG and "
                "the MPI-overlap Schwarz preconditioner");
    if(E->control.ala_geneo_preconditioner &&
       E->control.ala_two_level_preconditioner)
        myerror(E, "ALA GenEO and geometric two-level preconditioners are "
                "mutually exclusive");
    if(E->control.ala_feasibility_audit &&
       (!E->control.ala_pressure_buoyancy ||
        strcmp(E->control.uzawa,"ala_cg") != 0))
        myerror(E, "ala_feasibility_audit requires "
                "compressible_formulation=ala and uzawa=ala_cg");
    if(E->control.ala_element_vanka_smoother &&
       (strcmp(E->control.SOLVER_TYPE,"multigrid") != 0 ||
        !E->control.ala_pressure_buoyancy))
        myerror(E, "ala_element_vanka_smoother requires multigrid, "
                "and compressible_formulation=ala");

    PUTS(("\n"));

    Py_INCREF(Py_None);
    return Py_None;
}




/*==========================================================*/
/* helper functions */

FILE *get_output_stream(PyObject *out, struct All_variables*E)
{
    if (PyFile_Check(out)) {
        return PyFile_AsFile(out);
    }
    return NULL;
}


int _getStringProperty(PyObject* properties, char* attribute,
                       char* value, size_t valueSize, FILE* fp)
{
    PyObject *prop;
    char *buffer;
    Py_ssize_t length;

    if (!(prop = PyObject_GetAttrString(properties, attribute)))
        return -1;
    if (-1 == PyString_AsStringAndSize(prop, &buffer, &length))
        return -1;

    if (length >= (Py_ssize_t)valueSize) {
        PyErr_Format(PyExc_ValueError,
                     "value of '%s' cannot exceed %zu characters in length",
                     attribute, valueSize);
        return -1;
    }
    strcpy(value, buffer);

    if (fp)
        fprintf(fp, "%s=%s\n", attribute, value);

    return 0;
}


#define getTYPEProperty _getIntProperty
#define getTYPEVectorProperty _getIntVectorProperty
#define PyTYPE_Check PyInt_Check
#define CTYPE int
#define PyTYPE_AsCTYPE PyInt_AsLong
#define MESSAGE "an integer is required"
#define FORMAT "%d"
#include "getProperty.h"

#undef getTYPEProperty
#undef getTYPEVectorProperty
#undef PyTYPE_Check
#undef CTYPE
#undef PyTYPE_AsCTYPE
#undef MESSAGE
#undef FORMAT

#define getTYPEProperty _getFloatProperty
#define getTYPEVectorProperty _getFloatVectorProperty
#define PyTYPE_Check PyFloat_Check
#define CTYPE float
#define PyTYPE_AsCTYPE PyFloat_AsDouble
#define MESSAGE "a float is required"
#define FORMAT "%g"
#include "getProperty.h"


#undef getTYPEProperty
#undef getTYPEVectorProperty
#undef PyTYPE_Check
#undef CTYPE
#undef PyTYPE_AsCTYPE
#undef MESSAGE
#undef FORMAT

#define getTYPEProperty _getDoubleProperty
#define getTYPEVectorProperty _getDoubleVectorProperty
#define PyTYPE_Check PyFloat_Check
#define CTYPE double
#define PyTYPE_AsCTYPE PyFloat_AsDouble
#define MESSAGE "a float is required"
#define FORMAT "%g"
#include "getProperty.h"


/* $Id: setProperties.c 8157 2007-10-19 19:24:45Z tan2 $ */

/* End of file */
