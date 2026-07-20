/* Shared NODE/GP/ELEMENT depth-frequency profile engine. */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "npz_writer.h"
#include "output.h"
#include "parsing.h"

#define PROFILE_PERCENTILE_COUNT 5
#define PROFILE_STATS_COUNT 3
#define PROFILE_TEMPERATURE_BINS 200
#define PROFILE_TEMPERATURE_MIN_K 0.0
#define PROFILE_TEMPERATURE_MAX_K 4500.0

enum Profile_location {
    NODE_PROFILE = 0,
    GP_PROFILE = 1,
    ELEMENT_PROFILE = 2
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
    PROFILE_FLAG_TEMPERATURE = 0
};

extern void parallel_process_termination();
extern void get_global_shape_fn(struct All_variables *, int,
                                struct Shape_function *,
                                struct Shape_function_dx *,
                                struct Shape_function_dA *, int, int,
                                double [4][9], int, int);

static const double profile_percentile_levels[PROFILE_PERCENTILE_COUNT] =
    {10.0, 25.0, 50.0, 75.0, 90.0};

static double temperature_at_node(struct All_variables *E, int cap, int node)
{
    return E->data.Ttop + E->T[cap][node] * E->data.ref_temperature;
}

/* Add variables here. The statistical engine does not change. */
static const struct Profile_variable profile_variables[] = {
    {"temperature", "K", NODE_PROFILE, PROFILE_FLAG_TEMPERATURE,
     PROFILE_TEMPERATURE_BINS, PROFILE_TEMPERATURE_MIN_K,
     PROFILE_TEMPERATURE_MAX_K, temperature_at_node, NULL, NULL}
};

#define PROFILE_VARIABLE_COUNT \
    ((int)(sizeof(profile_variables) / sizeof(profile_variables[0])))

static int profile_flag(const struct All_variables *E, int flag_id)
{
    if (flag_id == PROFILE_FLAG_TEMPERATURE)
        return E->output.profile_temperature;
    return 0;
}

static void set_profile_flag(struct All_variables *E, int flag_id, int value)
{
    if (flag_id == PROFILE_FLAG_TEMPERATURE)
        E->output.profile_temperature = value;
}

static const struct Profile_variable *find_profile_variable(const char *name)
{
    int i;
    for (i = 0; i < PROFILE_VARIABLE_COUNT; ++i)
        if (strcmp(name, profile_variables[i].name) == 0)
            return &profile_variables[i];
    return NULL;
}

static double sample_profile_variable(struct All_variables *E,
                                      const struct Profile_variable *variable,
                                      int cap, int element, int gp)
{
    double value;
    int node_index, node;
    const int ends = enodes[E->mesh.nsd];

    if (variable->location == NODE_PROFILE) {
        value = 0.0;
        for (node_index = 1; node_index <= ends; ++node_index) {
            node = E->ien[cap][element].node[node_index];
            value += variable->node_getter(E, cap, node)
                * E->N.vpt[GNVINDEX(node_index, gp)];
        }
        return value;
    }
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

static void initialize_result(struct Profile_result *result,
                              const struct Profile_variable *variable,
                              int depth_bins, struct All_variables *E)
{
    int depth;
    const size_t histogram_size = (size_t)depth_bins * variable->bins;

    memset(result, 0, sizeof(*result));
    result->variable = variable;
    result->local_histogram = profile_allocate(histogram_size, sizeof(double), E);
    result->global_histogram = profile_allocate(histogram_size, sizeof(double), E);
    result->local_sum = profile_allocate(depth_bins, sizeof(double), E);
    result->global_sum = profile_allocate(depth_bins, sizeof(double), E);
    result->local_minimum = profile_allocate(depth_bins, sizeof(double), E);
    result->global_minimum = profile_allocate(depth_bins, sizeof(double), E);
    result->local_maximum = profile_allocate(depth_bins, sizeof(double), E);
    result->global_maximum = profile_allocate(depth_bins, sizeof(double), E);
    result->local_underflow = profile_allocate(depth_bins, sizeof(double), E);
    result->global_underflow = profile_allocate(depth_bins, sizeof(double), E);
    result->local_overflow = profile_allocate(depth_bins, sizeof(double), E);
    result->global_overflow = profile_allocate(depth_bins, sizeof(double), E);

    for (depth = 0; depth < depth_bins; ++depth) {
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

static void reduce_result(struct Profile_result *result, int depth_bins,
                          struct All_variables *E)
{
    const int histogram_size = depth_bins * result->variable->bins;

    MPI_Reduce(result->local_histogram, result->global_histogram,
               histogram_size, MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_sum, result->global_sum, depth_bins,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_minimum, result->global_minimum, depth_bins,
               MPI_DOUBLE, MPI_MIN, 0, E->parallel.world);
    MPI_Reduce(result->local_maximum, result->global_maximum, depth_bins,
               MPI_DOUBLE, MPI_MAX, 0, E->parallel.world);
    MPI_Reduce(result->local_underflow, result->global_underflow, depth_bins,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(result->local_overflow, result->global_overflow, depth_bins,
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
    double target, cumulative, within;
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
                return variable->minimum + (bin + 0.5) * width;
            within = (target - cumulative) / histogram[bin];
            return variable->minimum + (bin + within) * width;
        }
        cumulative += histogram[bin];
    }
    return variable->maximum;
}

static int add_f64(struct Npz_writer *writer, const char *name,
                   const double *data, int ndim, const size_t *shape)
{
    return npz_add_f64(writer, name, data, ndim, shape);
}

static int write_result(struct Npz_writer *writer,
                        const struct Profile_result *result,
                        const double *layer_volume, int depth_bins)
{
    const struct Profile_variable *variable = result->variable;
    double *bins, *histogram_fraction, *statistics, *percentiles;
    double *underflow_fraction, *overflow_fraction;
    char key[160];
    size_t shape_1d[1], shape_histogram[2], shape_stats[2];
    int depth, bin, p, status;

    bins = (double *)malloc((variable->bins + 1) * sizeof(double));
    histogram_fraction = (double *)calloc(
        (size_t)depth_bins * variable->bins, sizeof(double));
    statistics = (double *)calloc(
        (size_t)depth_bins * PROFILE_STATS_COUNT, sizeof(double));
    percentiles = (double *)calloc(
        (size_t)depth_bins * PROFILE_PERCENTILE_COUNT, sizeof(double));
    underflow_fraction = (double *)calloc(depth_bins, sizeof(double));
    overflow_fraction = (double *)calloc(depth_bins, sizeof(double));
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

    for (depth = 0; depth < depth_bins; ++depth) {
        if (layer_volume[depth] > 0.0) {
            for (bin = 0; bin < variable->bins; ++bin)
                histogram_fraction[depth * variable->bins + bin] =
                    result->global_histogram[depth * variable->bins + bin]
                    / layer_volume[depth];
            statistics[depth * PROFILE_STATS_COUNT] =
                result->global_minimum[depth];
            statistics[depth * PROFILE_STATS_COUNT + 1] =
                result->global_sum[depth] / layer_volume[depth];
            statistics[depth * PROFILE_STATS_COUNT + 2] =
                result->global_maximum[depth];
            underflow_fraction[depth] = result->global_underflow[depth]
                / layer_volume[depth];
            overflow_fraction[depth] = result->global_overflow[depth]
                / layer_volume[depth];
        }
        for (p = 0; p < PROFILE_PERCENTILE_COUNT; ++p)
            percentiles[depth * PROFILE_PERCENTILE_COUNT + p] =
                histogram_percentile(result, depth, layer_volume[depth],
                                     profile_percentile_levels[p]);
    }

    status = 0;
    shape_1d[0] = variable->bins + 1;
    snprintf(key, sizeof(key), "%s_bins", variable->name);
    status |= add_f64(writer, key, bins, 1, shape_1d);

    shape_histogram[0] = depth_bins;
    shape_histogram[1] = variable->bins;
    snprintf(key, sizeof(key), "%s_histogram_volume_km3", variable->name);
    status |= add_f64(writer, key, result->global_histogram, 2,
                      shape_histogram);
    snprintf(key, sizeof(key), "%s_hist2d", variable->name);
    status |= add_f64(writer, key, histogram_fraction, 2, shape_histogram);

    shape_stats[0] = depth_bins;
    shape_stats[1] = PROFILE_STATS_COUNT;
    snprintf(key, sizeof(key), "%s_stats", variable->name);
    status |= add_f64(writer, key, statistics, 2, shape_stats);

    shape_stats[1] = PROFILE_PERCENTILE_COUNT;
    snprintf(key, sizeof(key), "%s_percentiles", variable->name);
    status |= add_f64(writer, key, percentiles, 2, shape_stats);

    shape_1d[0] = depth_bins;
    snprintf(key, sizeof(key), "%s_underflow_fraction", variable->name);
    status |= add_f64(writer, key, underflow_fraction, 1, shape_1d);
    snprintf(key, sizeof(key), "%s_overflow_fraction", variable->name);
    status |= add_f64(writer, key, overflow_fraction, 1, shape_1d);

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
                              const double *layer_volume, int depth_bins)
{
    struct Npz_writer writer;
    char filename[255];
    double *depth_edges, *depth_mid;
    double elapsed_time, delta_temperature, surface_temperature;
    size_t shape[1];
    int schema_version, location_codes[PROFILE_VARIABLE_COUNT];
    int i, depth, global_z, status;

    depth_edges = (double *)malloc((depth_bins + 1) * sizeof(double));
    depth_mid = (double *)malloc(depth_bins * sizeof(double));
    if (!depth_edges || !depth_mid) {
        free(depth_edges);
        free(depth_mid);
        return -1;
    }

    for (depth = 0; depth <= depth_bins; ++depth) {
        global_z = depth_bins - depth;
        depth_edges[depth] = (E->sphere.ro - global_radii[global_z])
            * E->data.radius_km;
    }
    for (depth = 0; depth < depth_bins; ++depth)
        depth_mid[depth] = 0.5 * (depth_edges[depth] + depth_edges[depth + 1]);

    snprintf(filename, sizeof(filename), "%s/profiles_hist_%d.npz",
             E->control.data_dir, cycles);
    if (npz_open(&writer, filename) != 0) {
        fprintf(stderr, "Cannot open profile output '%s'\n", filename);
        free(depth_edges);
        free(depth_mid);
        return -1;
    }

    status = 0;
    schema_version = 1;
    status |= npz_add_i32(&writer, "schema_version", &schema_version, 0, NULL);
    status |= npz_add_i32(&writer, "step", &cycles, 0, NULL);
    elapsed_time = E->monitor.elapsed_time;
    delta_temperature = E->data.ref_temperature;
    surface_temperature = E->data.Ttop;
    status |= add_f64(&writer, "elapsed_time", &elapsed_time, 0, NULL);
    status |= add_f64(&writer, "DeltaT", &delta_temperature, 0, NULL);
    status |= add_f64(&writer, "T_surface_K", &surface_temperature, 0, NULL);

    shape[0] = depth_bins + 1;
    status |= add_f64(&writer, "depth_edges_km", depth_edges, 1, shape);
    shape[0] = depth_bins;
    status |= add_f64(&writer, "depth_km", depth_mid, 1, shape);
    status |= add_f64(&writer, "depth_volume_km3", layer_volume, 1, shape);
    shape[0] = PROFILE_PERCENTILE_COUNT;
    status |= add_f64(&writer, "percentile_levels",
                      profile_percentile_levels, 1, shape);

    for (i = 0; i < result_count; ++i) {
        location_codes[i] = (int)results[i].variable->location;
        snprintf(filename, sizeof(filename), "%s_location", results[i].variable->name);
        status |= npz_add_i32(&writer, filename, &location_codes[i], 0, NULL);
        status |= write_result(&writer, &results[i], layer_volume, depth_bins);
    }

    if (npz_close(&writer) != 0)
        status = -1;
    free(depth_edges);
    free(depth_mid);
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
    double *local_layer_volume, *global_layer_volume;
    double *local_radii, *global_radii;
    double value, volume, volume_scale;
    int enabled_indices[PROFILE_VARIABLE_COUNT];
    int result_count, variable_index, result_index;
    int depth_bins, cap, element, gp, depth, local_z, global_z;
    int radial_node, node, status;
    const int vpts = vpoints[E->mesh.nsd];
    const int sphere_key = 1;
    const int lev = E->mesh.levmax;

    result_count = 0;
    for (variable_index = 0; variable_index < PROFILE_VARIABLE_COUNT;
         ++variable_index)
        if (profile_flag(E, profile_variables[variable_index].flag_id))
            enabled_indices[result_count++] = variable_index;
    if (result_count == 0)
        return;

    depth_bins = E->mesh.noz - 1;
    volume_scale = E->data.radius_km * E->data.radius_km
        * E->data.radius_km;
    local_layer_volume = profile_allocate(depth_bins, sizeof(double), E);
    global_layer_volume = profile_allocate(depth_bins, sizeof(double), E);
    local_radii = profile_allocate(E->mesh.noz, sizeof(double), E);
    global_radii = profile_allocate(E->mesh.noz, sizeof(double), E);

    for (result_index = 0; result_index < result_count; ++result_index)
        initialize_result(&results[result_index],
                          &profile_variables[enabled_indices[result_index]],
                          depth_bins, E);

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
            depth = depth_bins - 1 - global_z;
            get_global_shape_fn(E, element, &GN, &GNx, &dOmega, 0,
                                sphere_key, rtf, lev, cap);

            for (gp = 1; gp <= vpts; ++gp) {
                volume = fabs(dOmega.vpt[gp]) * volume_scale;
                local_layer_volume[depth] += volume;
                for (result_index = 0; result_index < result_count;
                     ++result_index) {
                    value = sample_profile_variable(
                        E, results[result_index].variable, cap, element, gp);
                    accumulate_sample(&results[result_index], depth,
                                      value, volume);
                }
            }
        }

    MPI_Reduce(local_layer_volume, global_layer_volume, depth_bins,
               MPI_DOUBLE, MPI_SUM, 0, E->parallel.world);
    MPI_Reduce(local_radii, global_radii, E->mesh.noz, MPI_DOUBLE,
               MPI_MAX, 0, E->parallel.world);
    for (result_index = 0; result_index < result_count; ++result_index)
        reduce_result(&results[result_index], depth_bins, E);

    status = 0;
    if (E->parallel.me == 0)
        status = write_profiles_npz(E, cycles, results, result_count,
                                    global_radii, global_layer_volume,
                                    depth_bins);
    MPI_Bcast(&status, 1, MPI_INT, 0, E->parallel.world);
    if (status != 0) {
        if (E->parallel.me == 0)
            fprintf(stderr, "Failed to write profiles_hist_%d.npz\n", cycles);
        parallel_process_termination();
    }

    for (result_index = 0; result_index < result_count; ++result_index)
        free_result(&results[result_index]);
    free(local_layer_volume);
    free(global_layer_volume);
    free(local_radii);
    free(global_radii);
}
