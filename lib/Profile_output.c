/* Shared NODE/GP/ELEMENT depth-frequency profile engine. */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "material_properties.h"
#include "npz_writer.h"
#include "output.h"
#include "parsing.h"

#define PROFILE_PERCENTILE_COUNT 5
#define PROFILE_STATS_COUNT 3
#define PROFILE_TEMPERATURE_BINS 200
#define PROFILE_TEMPERATURE_MIN_K 0.0
#define PROFILE_TEMPERATURE_MAX_K 4500.0
#define PROFILE_K_BINS 200
#define PROFILE_K_MIN 0.0
#define PROFILE_K_MAX 35.0
#define PROFILE_VISCOSITY_BINS 200
#define PROFILE_VISCOSITY_MIN 18.0
#define PROFILE_VISCOSITY_MAX 24.0
#define PROFILE_DEFAULT_BINS 200
#define PROFILE_ALA_LOGABS_MIN -16.0
#define PROFILE_ALA_LOGABS_MAX 4.0
#define PROFILE_RELATIVE_MIN 0.0
#define PROFILE_RELATIVE_MAX 1.0000001
#define PROFILE_PRESSURE_MIN -3.0e9
#define PROFILE_PRESSURE_MAX 3.0e9
#define PROFILE_VELOCITY_MIN 0.0
#define PROFILE_VELOCITY_MAX 3.0e5
#define PROFILE_HEATING_MIN -1.0e4
#define PROFILE_HEATING_MAX 1.0e4
#define PROFILE_HEATING_POSITIVE_MAX 1.0e4
#define PROFILE_INTERNAL_MIN 0.0
#define PROFILE_INTERNAL_MAX 5.0e2
#define PROFILE_ASSIM_MIN -1.0e8
#define PROFILE_ASSIM_MAX 1.0e8
#define PROFILE_BUOYANCY_THERMAL_MIN -2.0e9
#define PROFILE_BUOYANCY_THERMAL_MAX 2.0e9
#define PROFILE_BUOYANCY_CHEMICAL_MIN -2.5e9
#define PROFILE_BUOYANCY_CHEMICAL_MAX 2.5e9
#define PROFILE_BUOYANCY_410_MIN -7.0e8
#define PROFILE_BUOYANCY_410_MAX 7.0e8
#define PROFILE_BUOYANCY_520_MIN -2.5e8
#define PROFILE_BUOYANCY_520_MAX 2.5e8
#define PROFILE_BUOYANCY_660_MIN -1.8e9
#define PROFILE_BUOYANCY_660_MAX 1.8e9
#define PROFILE_BUOYANCY_TOTAL_MIN -5.0e9
#define PROFILE_BUOYANCY_TOTAL_MAX 5.0e9

enum Profile_location {
    NODE_PROFILE = 0,    /* native node value, lumped nodal-volume weight */
    GP_PROFILE = 1,      /* native Gauss-point value */
    ELEMENT_PROFILE = 2 /* element value integrated over its Gauss points */
};

struct Profile_variable;

typedef double (*Profile_node_getter)(struct All_variables *, int, int);
typedef double (*Profile_gp_getter)(struct All_variables *, int, int, int);
typedef double (*Profile_element_getter)(struct All_variables *, int, int);

struct Profile_variable {
    const char *name;
    const char *units;
    enum Profile_location location;
    int flag_id;
    int bins;
    double minimum;
    double maximum;
    Profile_node_getter node_getter;
    Profile_gp_getter gp_getter;
    Profile_element_getter element_getter;
};

struct Profile_result {
    const struct Profile_variable *variable;
    int depth_count;
    double *local_histogram;
    double *global_histogram;
    double *local_sum;
    double *global_sum;
    double *local_minimum;
    double *global_minimum;
    double *local_maximum;
    double *global_maximum;
    double *local_underflow;
    double *global_underflow;
    double *local_overflow;
    double *global_overflow;
};

enum Profile_flag {
    PROFILE_FLAG_TEMPERATURE = 0,
    PROFILE_FLAG_K,
    PROFILE_FLAG_VISCOSITY_GP,
    PROFILE_FLAG_ALA_RESIDUAL_LOG10_ABS,
    PROFILE_FLAG_ALA_RELATIVE_RESIDUAL,
    PROFILE_FLAG_PRESSURE,
    PROFILE_FLAG_VELOCITY_MAGNITUDE,
    PROFILE_FLAG_QVISC,
    PROFILE_FLAG_QADI,
    PROFILE_FLAG_QADI_BASE,
    PROFILE_FLAG_QPHASE_ADI,
    PROFILE_FLAG_QINTERNAL,
    PROFILE_FLAG_QASSIM,
    PROFILE_FLAG_QTOTAL,
    PROFILE_FLAG_QVISC_MINUS_QADI_BASE,
    PROFILE_FLAG_BUOYANCY_THERMAL,
    PROFILE_FLAG_BUOYANCY_CHEMICAL,
    PROFILE_FLAG_BUOYANCY_PHASE_410,
    PROFILE_FLAG_BUOYANCY_PHASE_520,
    PROFILE_FLAG_BUOYANCY_PHASE_660,
    PROFILE_FLAG_BUOYANCY_TOTAL,
    PROFILE_FLAG_COUNT
};

static double **ala_residual_log10_abs_cache;
static double **ala_relative_residual_cache;
static double **buoyancy_thermal_cache;
static double **buoyancy_chemical_cache;
static double **buoyancy_phase_cache[PHASE_TRANSITIONS];
static double **buoyancy_total_cache;

extern void parallel_process_termination();
extern void get_global_shape_fn(struct All_variables *, int,
                                struct Shape_function *,
                                struct Shape_function_dx *,
                                struct Shape_function_dA *, int, int,
                                double [4][9], int, int);
extern void assemble_div_u(struct All_variables *, double **, double **, int);
extern void assemble_c_u(struct All_variables *, double **, double **, int);
extern void remove_horiz_ave2(struct All_variables *, double **);

static const double profile_percentile_levels[PROFILE_PERCENTILE_COUNT] =
    {10.0, 25.0, 50.0, 75.0, 90.0};

static double temperature_at_node(struct All_variables *E, int cap, int node)
{
    return E->data.Ttop + E->T[cap][node] * E->data.ref_temperature;
}

static double conductivity_at_gp(struct All_variables *E, int cap,
                                 int element, int gp)
{
    double temperature, kT, k_nd;
    int node_index, node;
    const int ends = enodes[E->mesh.nsd];

    temperature = 0.0;
    for (node_index = 1; node_index <= ends; ++node_index) {
        node = E->ien[cap][element].node[node_index];
        temperature += E->T[cap][node]
            * E->N.vpt[GNVINDEX(node_index, gp)];
    }
    kT = conductivity_temperature_factor(E, temperature);
    k_nd = conductivity_element_prefactor(
        E, cap, element, E->control.reference_conductivity) * kT;
    return k_nd * E->data.k0;
}

static double viscosity_at_gp(struct All_variables *E, int cap,
                              int element, int gp)
{
    double viscosity;
    const int vpts = vpoints[E->mesh.nsd];

    viscosity = E->EVI[E->mesh.levmax][cap][(element - 1) * vpts + gp]
        * E->data.ref_viscosity;
    if (viscosity <= 0.0 || !isfinite(viscosity)) {
        fprintf(stderr,
                "Invalid GP viscosity: cap=%d element=%d gp=%d eta=%e\n",
                cap, element, gp, viscosity);
        parallel_process_termination();
    }
    return log10(viscosity);
}

static double velocity_magnitude_at_node(struct All_variables *E,
                                         int cap, int node)
{
    double v1 = E->sphere.cap[cap].V[1][node];
    double v2 = E->sphere.cap[cap].V[2][node];
    double v3 = E->sphere.cap[cap].V[3][node];
    return sqrt(v1*v1 + v2*v2 + v3*v3);
}

static double pressure_at_element(struct All_variables *E, int cap,
                                  int element)
{
    return E->P[cap][element];
}

static double ala_residual_log10_abs_at_element(struct All_variables *E,
                                                int cap, int element)
{
    (void)E;
    return ala_residual_log10_abs_cache[cap][element];
}

static double ala_relative_residual_at_element(struct All_variables *E,
                                               int cap, int element)
{
    (void)E;
    return ala_relative_residual_cache[cap][element];
}

static double qvisc_at_element(struct All_variables *E, int cap, int element)
{
    return E->heating_visc[cap][element];
}

static double qadi_at_element(struct All_variables *E, int cap, int element)
{
    return E->heating_adi[cap][element];
}

static double qadi_base_at_element(struct All_variables *E, int cap,
                                   int element)
{
    return E->heating_adi_base[cap][element];
}

static double qphase_adi_at_element(struct All_variables *E, int cap,
                                    int element)
{
    return E->heating_phase_adi[cap][element];
}

static double qinternal_at_element(struct All_variables *E, int cap,
                                   int element)
{
    int radial_element = (element - 1) % E->lmesh.elz + 1;
    double rho = 0.5 * (E->refstate.rho[radial_element]
                        + E->refstate.rho[radial_element + 1]);
    double heat_production = E->control.Q0;

    if (E->control.tracer_enriched) {
        heat_production *= 1.0 - E->composition.comp_el[cap][0][element];
        heat_production += E->composition.comp_el[cap][0][element]
            * E->control.Q0ER;
    }
    return rho * heat_production;
}

static double qassim_at_element(struct All_variables *E, int cap, int element)
{
    return E->heating_assim[cap][element];
}

static double qtotal_at_element(struct All_variables *E, int cap, int element)
{
    return qinternal_at_element(E, cap, element)
        + E->heating_visc[cap][element] - E->heating_adi[cap][element]
        + E->heating_assim[cap][element];
}

static double qvisc_minus_qadi_base_at_element(struct All_variables *E,
                                               int cap, int element)
{
    return E->heating_visc[cap][element]
        - E->heating_adi_base[cap][element];
}

static double buoyancy_thermal_at_node(struct All_variables *E,
                                       int cap, int node)
{
    (void)E;
    return buoyancy_thermal_cache[cap][node];
}

static double buoyancy_chemical_at_node(struct All_variables *E,
                                        int cap, int node)
{
    (void)E;
    return buoyancy_chemical_cache[cap][node];
}

static double buoyancy_phase_410_at_node(struct All_variables *E,
                                         int cap, int node)
{
    (void)E;
    return buoyancy_phase_cache[0][cap][node];
}

static double buoyancy_phase_520_at_node(struct All_variables *E,
                                         int cap, int node)
{
    (void)E;
    return buoyancy_phase_cache[1][cap][node];
}

static double buoyancy_phase_660_at_node(struct All_variables *E,
                                         int cap, int node)
{
    (void)E;
    return buoyancy_phase_cache[2][cap][node];
}

static double buoyancy_total_at_node(struct All_variables *E,
                                     int cap, int node)
{
    (void)E;
    return buoyancy_total_cache[cap][node];
}

/* Add variables here. The statistical engine does not change. */
static const struct Profile_variable profile_variables[] = {
    {"temperature", "K", NODE_PROFILE, PROFILE_FLAG_TEMPERATURE,
     PROFILE_TEMPERATURE_BINS, PROFILE_TEMPERATURE_MIN_K,
     PROFILE_TEMPERATURE_MAX_K, temperature_at_node, NULL, NULL},
    {"k", "W/(m K)", GP_PROFILE, PROFILE_FLAG_K,
     PROFILE_K_BINS, PROFILE_K_MIN, PROFILE_K_MAX,
     NULL, conductivity_at_gp, NULL},
    {"viscosity_gp", "log10(Pa s)", GP_PROFILE,
     PROFILE_FLAG_VISCOSITY_GP, PROFILE_VISCOSITY_BINS,
     PROFILE_VISCOSITY_MIN, PROFILE_VISCOSITY_MAX,
     NULL, viscosity_at_gp, NULL},
    {"ala_residual_log10_abs", "log10(abs(Rala/V))", ELEMENT_PROFILE,
     PROFILE_FLAG_ALA_RESIDUAL_LOG10_ABS, PROFILE_DEFAULT_BINS,
     PROFILE_ALA_LOGABS_MIN, PROFILE_ALA_LOGABS_MAX,
     NULL, NULL, ala_residual_log10_abs_at_element},
    {"ala_relative_residual", "1", ELEMENT_PROFILE,
     PROFILE_FLAG_ALA_RELATIVE_RESIDUAL, PROFILE_DEFAULT_BINS,
     PROFILE_RELATIVE_MIN, PROFILE_RELATIVE_MAX,
     NULL, NULL, ala_relative_residual_at_element},
    {"pressure", "nondimensional", ELEMENT_PROFILE,
     PROFILE_FLAG_PRESSURE, PROFILE_DEFAULT_BINS,
     PROFILE_PRESSURE_MIN, PROFILE_PRESSURE_MAX,
     NULL, NULL, pressure_at_element},
    {"velocity_magnitude", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_VELOCITY_MAGNITUDE, PROFILE_DEFAULT_BINS,
     PROFILE_VELOCITY_MIN, PROFILE_VELOCITY_MAX,
     velocity_magnitude_at_node, NULL, NULL},
    {"qvisc", "nondimensional", ELEMENT_PROFILE, PROFILE_FLAG_QVISC,
     PROFILE_DEFAULT_BINS, 0.0, PROFILE_HEATING_POSITIVE_MAX,
     NULL, NULL, qvisc_at_element},
    {"qadi", "nondimensional", ELEMENT_PROFILE, PROFILE_FLAG_QADI,
     PROFILE_DEFAULT_BINS, PROFILE_HEATING_MIN, PROFILE_HEATING_MAX,
     NULL, NULL, qadi_at_element},
    {"qadi_base", "nondimensional", ELEMENT_PROFILE,
     PROFILE_FLAG_QADI_BASE, PROFILE_DEFAULT_BINS,
     PROFILE_HEATING_MIN, PROFILE_HEATING_MAX,
     NULL, NULL, qadi_base_at_element},
    {"qphase_adi", "nondimensional", ELEMENT_PROFILE,
     PROFILE_FLAG_QPHASE_ADI, PROFILE_DEFAULT_BINS,
     PROFILE_HEATING_MIN, PROFILE_HEATING_MAX,
     NULL, NULL, qphase_adi_at_element},
    {"qinternal", "nondimensional", ELEMENT_PROFILE,
     PROFILE_FLAG_QINTERNAL, PROFILE_DEFAULT_BINS,
     PROFILE_INTERNAL_MIN, PROFILE_INTERNAL_MAX,
     NULL, NULL, qinternal_at_element},
    {"qassim", "nondimensional", ELEMENT_PROFILE, PROFILE_FLAG_QASSIM,
     PROFILE_DEFAULT_BINS, PROFILE_ASSIM_MIN, PROFILE_ASSIM_MAX,
     NULL, NULL, qassim_at_element},
    {"qtotal", "nondimensional", ELEMENT_PROFILE, PROFILE_FLAG_QTOTAL,
     PROFILE_DEFAULT_BINS, PROFILE_ASSIM_MIN, PROFILE_ASSIM_MAX,
     NULL, NULL, qtotal_at_element},
    {"qvisc_minus_qadi_base", "nondimensional", ELEMENT_PROFILE,
     PROFILE_FLAG_QVISC_MINUS_QADI_BASE, PROFILE_DEFAULT_BINS,
     PROFILE_HEATING_MIN, PROFILE_HEATING_MAX,
     NULL, NULL, qvisc_minus_qadi_base_at_element},
    {"buoyancy_thermal", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_THERMAL, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_THERMAL_MIN, PROFILE_BUOYANCY_THERMAL_MAX,
     buoyancy_thermal_at_node, NULL, NULL},
    {"buoyancy_chemical", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_CHEMICAL, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_CHEMICAL_MIN, PROFILE_BUOYANCY_CHEMICAL_MAX,
     buoyancy_chemical_at_node, NULL, NULL},
    {"buoyancy_phase_410", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_PHASE_410, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_410_MIN, PROFILE_BUOYANCY_410_MAX,
     buoyancy_phase_410_at_node, NULL, NULL},
    {"buoyancy_phase_520", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_PHASE_520, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_520_MIN, PROFILE_BUOYANCY_520_MAX,
     buoyancy_phase_520_at_node, NULL, NULL},
    {"buoyancy_phase_660", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_PHASE_660, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_660_MIN, PROFILE_BUOYANCY_660_MAX,
     buoyancy_phase_660_at_node, NULL, NULL},
    {"buoyancy_total", "nondimensional", NODE_PROFILE,
     PROFILE_FLAG_BUOYANCY_TOTAL, PROFILE_DEFAULT_BINS,
     PROFILE_BUOYANCY_TOTAL_MIN, PROFILE_BUOYANCY_TOTAL_MAX,
     buoyancy_total_at_node, NULL, NULL}
};

#define PROFILE_VARIABLE_COUNT \
    ((int)(sizeof(profile_variables) / sizeof(profile_variables[0])))

static int profile_flag(const struct All_variables *E, int flag_id)
{
    return flag_id >= 0 && flag_id < PROFILE_FLAG_COUNT
        ? E->output.profile_flags[flag_id] : 0;
}

static void set_profile_flag(struct All_variables *E, int flag_id, int value)
{
    if (flag_id >= 0 && flag_id < PROFILE_FLAG_COUNT)
        E->output.profile_flags[flag_id] = value;
}

static const struct Profile_variable *find_profile_variable(const char *name)
{
    int i;
    for (i = 0; i < PROFILE_VARIABLE_COUNT; ++i)
        if (strcmp(name, profile_variables[i].name) == 0)
            return &profile_variables[i];
    return NULL;
}

static double sample_gp_or_element(struct All_variables *E,
                                   const struct Profile_variable *variable,
                                   int cap, int element, int gp)
{
    if (variable->location == GP_PROFILE)
        return variable->gp_getter(E, cap, element, gp);
    if (variable->location == ELEMENT_PROFILE)
        return variable->element_getter(E, cap, element);

    return 0.0;
}

static void *profile_allocate(size_t count, size_t size,
                              struct All_variables *E)
{
    void *memory = calloc(count, size);
    if (!memory) {
        fprintf(stderr, "rank %d: cannot allocate profile arrays\n",
                E->parallel.me);
        parallel_process_termination();
    }
    return memory;
}

static double **allocate_profile_field(struct All_variables *E, int count)
{
    int cap;
    double **field = profile_allocate(NCS, sizeof(double *), E);

    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        field[cap] = profile_allocate(count, sizeof(double), E);
    return field;
}

static void free_profile_field(struct All_variables *E, double **field)
{
    int cap;

    if (!field)
        return;
    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        free(field[cap]);
    free(field);
}

static int any_profile_flag(struct All_variables *E, int first, int last)
{
    int flag;
    for (flag = first; flag <= last; ++flag)
        if (profile_flag(E, flag))
            return 1;
    return 0;
}

static void prepare_ala_profile_fields(struct All_variables *E)
{
    double **rdiv, **rbeta;
    double volume, residual, denominator, relative;
    int cap, element;

    if (!any_profile_flag(E, PROFILE_FLAG_ALA_RESIDUAL_LOG10_ABS,
                          PROFILE_FLAG_ALA_RELATIVE_RESIDUAL))
        return;

    ala_residual_log10_abs_cache = allocate_profile_field(
        E, E->lmesh.nel + 1);
    ala_relative_residual_cache = allocate_profile_field(
        E, E->lmesh.nel + 1);
    rdiv = allocate_profile_field(E, E->lmesh.nel + 1);
    rbeta = allocate_profile_field(E, E->lmesh.nel + 1);

    assemble_div_u(E, E->U, rdiv, E->mesh.levmax);
    assemble_c_u(E, E->U, rbeta, E->mesh.levmax);
    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        for (element = 1; element <= E->lmesh.nel; ++element) {
            volume = E->eco[cap][element].area;
            residual = volume > 0.0
                ? (rdiv[cap][element] + rbeta[cap][element]) / volume
                : 0.0;
            denominator = fabs(rdiv[cap][element])
                + fabs(rbeta[cap][element]);
            relative = fabs(rdiv[cap][element] + rbeta[cap][element])
                / max(denominator, 1.0e-30 * max(volume, 1.0));
            ala_residual_log10_abs_cache[cap][element] =
                log10(max(fabs(residual), 1.0e-16));
            ala_relative_residual_cache[cap][element] = min(relative, 1.0);
        }

    free_profile_field(E, rdiv);
    free_profile_field(E, rbeta);
}

static void prepare_buoyancy_profile_fields(struct All_variables *E)
{
    int cap, node, component, phase_index, radial_node;
    double gravity;

    if (!any_profile_flag(E, PROFILE_FLAG_BUOYANCY_THERMAL,
                          PROFILE_FLAG_BUOYANCY_TOTAL))
        return;

    buoyancy_thermal_cache = allocate_profile_field(E, E->lmesh.nno + 1);
    buoyancy_chemical_cache = allocate_profile_field(E, E->lmesh.nno + 1);
    buoyancy_total_cache = allocate_profile_field(E, E->lmesh.nno + 1);
    for (phase_index = 0; phase_index < PHASE_TRANSITIONS; ++phase_index)
        buoyancy_phase_cache[phase_index] = allocate_profile_field(
            E, E->lmesh.nno + 1);

    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        for (node = 1; node <= E->lmesh.nno; ++node) {
            radial_node = (node - 1) % E->lmesh.noz + 1;
            buoyancy_thermal_cache[cap][node] = E->control.Atemp
                * E->refstate.rho[radial_node]
                * E->refstate.thermal_expansivity[radial_node]
                * E->T[cap][node];
            if (E->control.tracer && E->composition.ichemical_buoyancy)
                for (component = 0; component < E->composition.ncomp;
                     ++component)
                    buoyancy_chemical_cache[cap][node] -= E->control.Atemp
                        * E->composition.buoyancy_ratio[component]
                        * E->composition.comp_node[cap][component][node];
            for (phase_index = 0; phase_index < PHASE_TRANSITIONS;
                 ++phase_index)
                if (E->control.phase[phase_index].Ra != 0.0)
                    buoyancy_phase_cache[phase_index][cap][node] =
                        -E->control.phase[phase_index].Ra
                        * E->phase_B[phase_index][cap][node];
        }

    remove_horiz_ave2(E, buoyancy_thermal_cache);
    remove_horiz_ave2(E, buoyancy_chemical_cache);
    for (phase_index = 0; phase_index < PHASE_TRANSITIONS; ++phase_index)
        remove_horiz_ave2(E, buoyancy_phase_cache[phase_index]);

    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        for (node = 1; node <= E->lmesh.nno; ++node) {
            radial_node = (node - 1) % E->lmesh.noz + 1;
            gravity = E->refstate.gravity[radial_node];
            buoyancy_thermal_cache[cap][node] *= gravity;
            buoyancy_chemical_cache[cap][node] *= gravity;
            buoyancy_total_cache[cap][node] =
                buoyancy_thermal_cache[cap][node]
                + buoyancy_chemical_cache[cap][node];
            for (phase_index = 0; phase_index < PHASE_TRANSITIONS;
                 ++phase_index) {
                buoyancy_phase_cache[phase_index][cap][node] *= gravity;
                buoyancy_total_cache[cap][node] +=
                    buoyancy_phase_cache[phase_index][cap][node];
            }
        }
}

static void prepare_derived_profile_fields(struct All_variables *E)
{
    prepare_ala_profile_fields(E);
    prepare_buoyancy_profile_fields(E);
}

static void free_derived_profile_fields(struct All_variables *E)
{
    int phase_index;

    free_profile_field(E, ala_residual_log10_abs_cache);
    free_profile_field(E, ala_relative_residual_cache);
    free_profile_field(E, buoyancy_thermal_cache);
    free_profile_field(E, buoyancy_chemical_cache);
    free_profile_field(E, buoyancy_total_cache);
    for (phase_index = 0; phase_index < PHASE_TRANSITIONS; ++phase_index)
        free_profile_field(E, buoyancy_phase_cache[phase_index]);
    ala_residual_log10_abs_cache = NULL;
    ala_relative_residual_cache = NULL;
    buoyancy_thermal_cache = NULL;
    buoyancy_chemical_cache = NULL;
    buoyancy_total_cache = NULL;
    for (phase_index = 0; phase_index < PHASE_TRANSITIONS; ++phase_index)
        buoyancy_phase_cache[phase_index] = NULL;
}

static void initialize_result(struct Profile_result *result,
                              const struct Profile_variable *variable,
                              struct All_variables *E)
{
    int depth;
    int depth_count;
    size_t histogram_size;

    memset(result, 0, sizeof(*result));
    result->variable = variable;
    depth_count = variable->location == NODE_PROFILE
        ? E->mesh.noz : E->mesh.noz - 1;
    result->depth_count = depth_count;
    histogram_size = (size_t)depth_count * variable->bins;
    result->local_histogram = profile_allocate(histogram_size, sizeof(double), E);
    result->global_histogram = profile_allocate(histogram_size, sizeof(double), E);
    result->local_sum = profile_allocate(depth_count, sizeof(double), E);
    result->global_sum = profile_allocate(depth_count, sizeof(double), E);
    result->local_minimum = profile_allocate(depth_count, sizeof(double), E);
    result->global_minimum = profile_allocate(depth_count, sizeof(double), E);
    result->local_maximum = profile_allocate(depth_count, sizeof(double), E);
    result->global_maximum = profile_allocate(depth_count, sizeof(double), E);
    result->local_underflow = profile_allocate(depth_count, sizeof(double), E);
    result->global_underflow = profile_allocate(depth_count, sizeof(double), E);
    result->local_overflow = profile_allocate(depth_count, sizeof(double), E);
    result->global_overflow = profile_allocate(depth_count, sizeof(double), E);

    for (depth = 0; depth < depth_count; ++depth) {
        result->local_minimum[depth] = DBL_MAX;
        result->local_maximum[depth] = -DBL_MAX;
    }
}

static void free_result(struct Profile_result *result)
{
    free(result->local_histogram);
    free(result->global_histogram);
    free(result->local_sum);
    free(result->global_sum);
    free(result->local_minimum);
    free(result->global_minimum);
    free(result->local_maximum);
    free(result->global_maximum);
    free(result->local_underflow);
    free(result->global_underflow);
    free(result->local_overflow);
    free(result->global_overflow);
}

static void accumulate_sample(struct Profile_result *result, int depth,
                              double value, double volume)
{
    const struct Profile_variable *variable = result->variable;
    const double width = (variable->maximum - variable->minimum)
        / variable->bins;
    int bin;

    result->local_sum[depth] += value * volume;
    if (value < result->local_minimum[depth])
        result->local_minimum[depth] = value;
    if (value > result->local_maximum[depth])
        result->local_maximum[depth] = value;

    if (value < variable->minimum)
        result->local_underflow[depth] += volume;
    else if (value >= variable->maximum)
        result->local_overflow[depth] += volume;
    else {
        bin = (int)((value - variable->minimum) / width);
        result->local_histogram[depth * variable->bins + bin] += volume;
    }
}

static void reduce_result(struct Profile_result *result,
                          struct All_variables *E)
{
    const int depth_count = result->depth_count;
    const int histogram_size = depth_count * result->variable->bins;

    MPI_Reduce(result->local_histogram, result->global_histogram,
               histogram_size, MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_sum, result->global_sum, depth_count,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_minimum, result->global_minimum, depth_count,
               MPI_DOUBLE, MPI_MIN, 0, E->parallel.world);
    MPI_Reduce(result->local_maximum, result->global_maximum, depth_count,
               MPI_DOUBLE, MPI_MAX, 0, E->parallel.world);
    MPI_Reduce(result->local_underflow, result->global_underflow, depth_count,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_overflow, result->global_overflow, depth_count,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
}

static double histogram_percentile(const struct Profile_result *result,
                                   int depth, double layer_volume,
                                   double percentile)
{
    const struct Profile_variable *variable = result->variable;
    const double *histogram = result->global_histogram
        + depth * variable->bins;
    const double width = (variable->maximum - variable->minimum)
        / variable->bins;
    double target, cumulative, within, value;
    int bin;

    if (layer_volume <= 0.0)
        return 0.0;
    target = percentile * 0.01 * layer_volume;
    if (target <= result->global_underflow[depth])
        return variable->minimum;

    cumulative = result->global_underflow[depth];
    for (bin = 0; bin < variable->bins; ++bin) {
        if (target <= cumulative + histogram[bin]) {
            if (histogram[bin] <= 0.0)
                value = variable->minimum + (bin + 0.5) * width;
            else {
                within = (target - cumulative) / histogram[bin];
                value = variable->minimum + (bin + within) * width;
            }
            if (value < result->global_minimum[depth])
                value = result->global_minimum[depth];
            if (value > result->global_maximum[depth])
                value = result->global_maximum[depth];
            return value;
        }
        cumulative += histogram[bin];
    }
    return result->global_maximum[depth];
}

static int add_f64(struct Npz_writer *writer, const char *name,
                   const double *data, int ndim, const size_t *shape)
{
    return npz_add_f64(writer, name, data, ndim, shape);
}

static int write_result(struct Npz_writer *writer,
                        const struct Profile_result *result,
                        const double *depth_coordinates,
                        const double *depth_volume)
{
    const struct Profile_variable *variable = result->variable;
    const int depth_count = result->depth_count;
    double *bins, *histogram_fraction, *statistics, *percentiles;
    double *underflow_fraction, *overflow_fraction;
    char key[160];
    size_t shape_1d[1], shape_histogram[2], shape_stats[2];
    int depth, bin, p, status;

    bins = (double *)malloc((variable->bins + 1) * sizeof(double));
    histogram_fraction = (double *)calloc(
        (size_t)depth_count * variable->bins, sizeof(double));
    statistics = (double *)calloc(
        (size_t)depth_count * PROFILE_STATS_COUNT, sizeof(double));
    percentiles = (double *)calloc(
        (size_t)depth_count * PROFILE_PERCENTILE_COUNT, sizeof(double));
    underflow_fraction = (double *)calloc(depth_count, sizeof(double));
    overflow_fraction = (double *)calloc(depth_count, sizeof(double));
    if (!bins || !histogram_fraction || !statistics || !percentiles ||
        !underflow_fraction || !overflow_fraction) {
        free(bins);
        free(histogram_fraction);
        free(statistics);
        free(percentiles);
        free(underflow_fraction);
        free(overflow_fraction);
        return -1;
    }

    for (bin = 0; bin <= variable->bins; ++bin)
        bins[bin] = variable->minimum
            + (variable->maximum - variable->minimum) * bin / variable->bins;

    for (depth = 0; depth < depth_count; ++depth) {
        if (depth_volume[depth] > 0.0) {
            for (bin = 0; bin < variable->bins; ++bin)
                histogram_fraction[depth * variable->bins + bin] =
                    result->global_histogram[depth * variable->bins + bin]
                    / depth_volume[depth];
            statistics[depth * PROFILE_STATS_COUNT] =
                result->global_minimum[depth];
            statistics[depth * PROFILE_STATS_COUNT + 1] =
                result->global_sum[depth] / depth_volume[depth];
            statistics[depth * PROFILE_STATS_COUNT + 2] =
                result->global_maximum[depth];
            underflow_fraction[depth] = result->global_underflow[depth]
                / depth_volume[depth];
            overflow_fraction[depth] = result->global_overflow[depth]
                / depth_volume[depth];
        }
        for (p = 0; p < PROFILE_PERCENTILE_COUNT; ++p)
            percentiles[depth * PROFILE_PERCENTILE_COUNT + p] =
                histogram_percentile(result, depth, depth_volume[depth],
                                     profile_percentile_levels[p]);
    }

    status = 0;
    shape_1d[0] = variable->bins + 1;
    snprintf(key, sizeof(key), "%s_bins", variable->name);
    status |= add_f64(writer, key, bins, 1, shape_1d);

    shape_histogram[0] = depth_count;
    shape_histogram[1] = variable->bins;
    snprintf(key, sizeof(key), "%s_histogram_volume_km3", variable->name);
    status |= add_f64(writer, key, result->global_histogram, 2,
                      shape_histogram);
    snprintf(key, sizeof(key), "%s_hist2d", variable->name);
    status |= add_f64(writer, key, histogram_fraction, 2, shape_histogram);

    shape_stats[0] = depth_count;
    shape_stats[1] = PROFILE_STATS_COUNT;
    snprintf(key, sizeof(key), "%s_stats", variable->name);
    status |= add_f64(writer, key, statistics, 2, shape_stats);

    shape_stats[1] = PROFILE_PERCENTILE_COUNT;
    snprintf(key, sizeof(key), "%s_percentiles", variable->name);
    status |= add_f64(writer, key, percentiles, 2, shape_stats);

    shape_1d[0] = depth_count;
    snprintf(key, sizeof(key), "%s_underflow_fraction", variable->name);
    status |= add_f64(writer, key, underflow_fraction, 1, shape_1d);
    snprintf(key, sizeof(key), "%s_overflow_fraction", variable->name);
    status |= add_f64(writer, key, overflow_fraction, 1, shape_1d);
    snprintf(key, sizeof(key), "%s_depth_km", variable->name);
    status |= add_f64(writer, key, depth_coordinates, 1, shape_1d);
    snprintf(key, sizeof(key), "%s_depth_volume_km3", variable->name);
    status |= add_f64(writer, key, depth_volume, 1, shape_1d);

    /* Compatibility aliases used by scripts/make_profiles_hist.py. */
    if (strcmp(variable->name, "temperature") == 0) {
        shape_1d[0] = variable->bins + 1;
        status |= add_f64(writer, "T_bins", bins, 1, shape_1d);
        status |= add_f64(writer, "hist2d", histogram_fraction, 2,
                          shape_histogram);
        shape_stats[1] = PROFILE_STATS_COUNT;
        status |= add_f64(writer, "stats", statistics, 2, shape_stats);
        shape_stats[1] = PROFILE_PERCENTILE_COUNT;
        status |= add_f64(writer, "percentiles", percentiles, 2,
                          shape_stats);
    }

    free(bins);
    free(histogram_fraction);
    free(statistics);
    free(percentiles);
    free(underflow_fraction);
    free(overflow_fraction);
    return status;
}

static int write_profiles_npz(struct All_variables *E, int cycles,
                              struct Profile_result *results,
                              int result_count, const double *global_radii,
                              const double *node_volume,
                              const double *element_volume)
{
    struct Npz_writer writer;
    char filename[255];
    double *node_depth, *element_depth;
    double elapsed_time, delta_temperature, surface_temperature;
    size_t shape[1];
    int schema_version, location_codes[PROFILE_VARIABLE_COUNT];
    int i, depth, global_z, status;
    const int node_depth_count = E->mesh.noz;
    const int element_depth_count = E->mesh.noz - 1;

    node_depth = (double *)malloc(node_depth_count * sizeof(double));
    element_depth = (double *)malloc(element_depth_count * sizeof(double));
    if (!node_depth || !element_depth) {
        free(node_depth);
        free(element_depth);
        return -1;
    }

    for (depth = 0; depth < node_depth_count; ++depth) {
        global_z = node_depth_count - 1 - depth;
        node_depth[depth] = (E->sphere.ro - global_radii[global_z])
            * E->data.radius_km;
    }
    for (depth = 0; depth < element_depth_count; ++depth)
        element_depth[depth] = 0.5
            * (node_depth[depth] + node_depth[depth + 1]);

    snprintf(filename, sizeof(filename), "%s/profiles_hist_%d.npz",
             E->control.data_dir, cycles);
    if (npz_open(&writer, filename) != 0) {
        fprintf(stderr, "Cannot open profile output '%s'\n", filename);
        free(node_depth);
        free(element_depth);
        return -1;
    }

    status = 0;
    schema_version = 2;
    status |= npz_add_i32(&writer, "schema_version", &schema_version, 0, NULL);
    status |= npz_add_i32(&writer, "step", &cycles, 0, NULL);
    elapsed_time = E->monitor.elapsed_time;
    delta_temperature = E->data.ref_temperature;
    surface_temperature = E->data.Ttop;
    status |= add_f64(&writer, "elapsed_time", &elapsed_time, 0, NULL);
    status |= add_f64(&writer, "DeltaT", &delta_temperature, 0, NULL);
    status |= add_f64(&writer, "T_surface_K", &surface_temperature, 0, NULL);

    shape[0] = node_depth_count;
    status |= add_f64(&writer, "node_depth_km", node_depth, 1, shape);
    status |= add_f64(&writer, "node_volume_km3", node_volume, 1, shape);
    status |= add_f64(&writer, "element_depth_edges_km", node_depth, 1, shape);
    shape[0] = element_depth_count;
    status |= add_f64(&writer, "element_depth_km", element_depth, 1, shape);
    status |= add_f64(&writer, "element_volume_km3", element_volume, 1, shape);

    /* Compatibility aliases follow the current NODE temperature profile. */
    shape[0] = node_depth_count;
    status |= add_f64(&writer, "depth_km", node_depth, 1, shape);
    status |= add_f64(&writer, "depth_volume_km3", node_volume, 1, shape);
    shape[0] = PROFILE_PERCENTILE_COUNT;
    status |= add_f64(&writer, "percentile_levels",
                      profile_percentile_levels, 1, shape);

    for (i = 0; i < result_count; ++i) {
        location_codes[i] = (int)results[i].variable->location;
        snprintf(filename, sizeof(filename), "%s_location", results[i].variable->name);
        status |= npz_add_i32(&writer, filename, &location_codes[i], 0, NULL);
        if (results[i].variable->location == NODE_PROFILE)
            status |= write_result(&writer, &results[i], node_depth,
                                   node_volume);
        else
            status |= write_result(&writer, &results[i], element_depth,
                                   element_volume);
    }

    if (npz_close(&writer) != 0)
        status = -1;
    free(node_depth);
    free(element_depth);
    return status;
}

void profile_input(struct All_variables *E)
{
    input_string("profile_optional", E->output.profile_optional, "",
                 E->parallel.me);
}

void profile_parse_optional(struct All_variables *E)
{
    char *strip(char *);
    const struct Profile_variable *variable;
    char *field, *next;
    int i;

    for (i = 0; i < PROFILE_VARIABLE_COUNT; ++i)
        set_profile_flag(E, profile_variables[i].flag_id, 0);

    next = E->output.profile_optional;
    while ((field = strsep(&next, ",")) != NULL) {
        field = strip(field);
        if (field[0] == '\0')
            continue;
        variable = find_profile_variable(field);
        if (variable)
            set_profile_flag(E, variable->flag_id, 1);
        else if (E->parallel.me == 0)
            fprintf(stderr,
                    "Warning: unknown field for profile_optional: %s\n",
                    field);
    }
}

void output_profiles(struct All_variables *E, int cycles)
{
    struct Profile_result results[PROFILE_VARIABLE_COUNT];
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    double rtf[4][9];
    double *local_node_volume, *global_node_volume;
    double *local_element_volume, *global_element_volume;
    double *local_radii, *global_radii;
    double value, volume, nodal_volume, volume_scale;
    int enabled_indices[PROFILE_VARIABLE_COUNT];
    int result_count, variable_index, result_index;
    int node_depth_count, element_depth_count;
    int cap, element, gp, element_depth, node_depth, local_z, global_z;
    int radial_node, local_radial_node, node_index, node, status;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int sphere_key = 1;
    const int lev = E->mesh.levmax;

    result_count = 0;
    for (variable_index = 0; variable_index < PROFILE_VARIABLE_COUNT;
         ++variable_index)
        if (profile_flag(E, profile_variables[variable_index].flag_id))
            enabled_indices[result_count++] = variable_index;
    if (result_count == 0)
        return;

    node_depth_count = E->mesh.noz;
    element_depth_count = E->mesh.noz - 1;
    volume_scale = E->data.radius_km * E->data.radius_km
        * E->data.radius_km;
    local_node_volume = profile_allocate(node_depth_count, sizeof(double), E);
    global_node_volume = profile_allocate(node_depth_count, sizeof(double), E);
    local_element_volume = profile_allocate(element_depth_count, sizeof(double), E);
    global_element_volume = profile_allocate(element_depth_count, sizeof(double), E);
    local_radii = profile_allocate(E->mesh.noz, sizeof(double), E);
    global_radii = profile_allocate(E->mesh.noz, sizeof(double), E);

    for (result_index = 0; result_index < result_count; ++result_index)
        initialize_result(&results[result_index],
                          &profile_variables[enabled_indices[result_index]],
                          E);
    prepare_derived_profile_fields(E);

    for (radial_node = 1; radial_node <= E->lmesh.noz; ++radial_node) {
        global_z = E->lmesh.nzs + radial_node - 2;
        node = radial_node;
        if (global_z >= 0 && global_z < E->mesh.noz)
            local_radii[global_z] = E->sx[1][3][node];
    }

    for (cap = 1; cap <= E->sphere.caps_per_proc; ++cap)
        for (element = 1; element <= E->lmesh.nel; ++element) {
            local_z = (element - 1) % E->lmesh.elz;
            global_z = E->lmesh.ezs + local_z;
            element_depth = element_depth_count - 1 - global_z;
            get_global_shape_fn(E, element, &GN, &GNx, &dOmega, 0,
                                sphere_key, rtf, lev, cap);

            for (gp = 1; gp <= vpts; ++gp) {
                volume = fabs(dOmega.vpt[gp]) * volume_scale;
                local_element_volume[element_depth] += volume;

                for (node_index = 1; node_index <= ends; ++node_index) {
                    node = E->ien[cap][element].node[node_index];
                    local_radial_node = (node - 1) % E->lmesh.noz + 1;
                    global_z = E->lmesh.nzs + local_radial_node - 2;
                    node_depth = node_depth_count - 1 - global_z;
                    nodal_volume = E->N.vpt[GNVINDEX(node_index, gp)]
                        * volume;
                    local_node_volume[node_depth] += nodal_volume;
                    for (result_index = 0; result_index < result_count;
                         ++result_index) {
                        if (results[result_index].variable->location
                            != NODE_PROFILE)
                            continue;
                        value = results[result_index].variable->node_getter(
                            E, cap, node);
                        accumulate_sample(&results[result_index], node_depth,
                                          value, nodal_volume);
                    }
                }

                for (result_index = 0; result_index < result_count;
                     ++result_index) {
                    if (results[result_index].variable->location
                        == NODE_PROFILE)
                        continue;
                    value = sample_gp_or_element(
                        E, results[result_index].variable, cap, element, gp);
                    accumulate_sample(&results[result_index], element_depth,
                                      value, volume);
                }
            }
        }

    MPI_Reduce(local_node_volume, global_node_volume, node_depth_count,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(local_element_volume, global_element_volume,
               element_depth_count,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(local_radii, global_radii, E->mesh.noz, MPI_DOUBLE,
               MPI_MAX, 0, E->parallel.world);
    for (result_index = 0; result_index < result_count; ++result_index)
        reduce_result(&results[result_index], E);

    status = 0;
    if (E->parallel.me == 0)
        status = write_profiles_npz(E, cycles, results, result_count,
                                    global_radii, global_node_volume,
                                    global_element_volume);
    MPI_Bcast(&status, 1, MPI_INT, 0, E->parallel.world);
    if (status != 0) {
        if (E->parallel.me == 0)
            fprintf(stderr, "Failed to write profiles_hist_%d.npz\n", cycles);
        parallel_process_termination();
    }

    for (result_index = 0; result_index < result_count; ++result_index)
        free_result(&results[result_index]);
    free_derived_profile_fields(E);
    free(local_node_volume);
    free(global_node_volume);
    free(local_element_volume);
    free(global_element_volume);
    free(local_radii);
    free(global_radii);
}
