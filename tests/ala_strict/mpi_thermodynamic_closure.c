/*
 * Standalone synthetic MPI thermal-projection smoke test.
 *
 * This test deliberately does not initialize CitcomS, use its mesh
 * quadrature, or call remove_horiz_ave2().  Its deterministic spherical-cell
 * weights exercise weighted MPI reduction algebra only; it is not a
 * production-MPI, phase, or composition validation.
 */

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_RADIAL_NODES 4096
#define SYNTHETIC_LATITUDE_CELLS 12
#define SYNTHETIC_LONGITUDE_CELLS 32
#define SYNTHETIC_SURFACE_CELLS \
    (SYNTHETIC_LATITUDE_CELLS * SYNTHETIC_LONGITUDE_CELLS)

static const double pi = 3.14159265358979323846264338327950288;
static const double coefficient_tolerance = 1.0e-10;
static const double scaled_relative_tolerance = 1.0e-10;
static const double minimum_signal = 1.0e-8;

struct norms {
    double weighted_l2;
    double weighted_rms;
    double linf;
    double total_weight;
    int finite;
};

static int load_reference(const char *path, double *rho, double *gravity,
                          double *tref, double *alpha)
{
    FILE *stream;
    char line[4096];
    int count = 0;

    stream = fopen(path, "r");
    if (stream == NULL)
        return -1;

    while (fgets(line, sizeof(line), stream) != NULL) {
        double values[8];
        char *cursor = line;
        char *end;
        int column;
        int columns;

        while (isspace((unsigned char)*cursor))
            ++cursor;
        if (*cursor == '\0' || *cursor == '#')
            continue;

        columns = 0;
        while (columns < 8) {
            errno = 0;
            values[columns] = strtod(cursor, &end);
            if (end == cursor)
                break;
            if (errno != 0) {
                fclose(stream);
                return -2;
            }
            cursor = end;
            ++columns;
        }
        while (isspace((unsigned char)*cursor))
            ++cursor;
        if (columns != 7 && columns != 8) {
            fclose(stream);
            return -2;
        }
        if (*cursor != '\0' && *cursor != '#') {
            fclose(stream);
            return -2;
        }
        if (count >= MAX_RADIAL_NODES) {
            fclose(stream);
            return -3;
        }
        for (column = 0; column < columns; ++column)
            if (!isfinite(values[column])) {
                fclose(stream);
                return -4;
            }
        if (values[0] <= 0.0 || values[1] <= 0.0
            || values[3] <= 0.0) {
            fclose(stream);
            return -4;
        }
        rho[count] = values[0];
        gravity[count] = values[1];
        tref[count] = values[2];
        alpha[count] = values[3];
        ++count;
    }
    fclose(stream);
    return count >= 2 ? count : -5;
}

static int parse_atemp(const char *text, double *requested, double *effective)
{
    char *end;
    double value;
    float production_storage;

    errno = 0;
    value = strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0'
        || !isfinite(value) || value <= 0.0)
        return 0;

    /* E->control.Atemp is float in production; preserve that rounding. */
    production_storage = (float)value;
    if (!isfinite(production_storage) || production_storage <= 0.0f)
        return 0;
    *requested = value;
    *effective = (double)production_storage;
    return 1;
}

static void partition_surface(int rank, int ranks, int *start, int *count)
{
    int quotient = SYNTHETIC_SURFACE_CELLS / ranks;
    int remainder = SYNTHETIC_SURFACE_CELLS % ranks;

    *count = quotient + (rank < remainder ? 1 : 0);
    *start = rank * quotient + (rank < remainder ? rank : remainder);
}

static double surface_weight(int global_cell)
{
    int latitude_cell = global_cell / SYNTHETIC_LONGITUDE_CELLS;
    double latitude_lower = -0.5 * pi
        + pi * (double)latitude_cell
        / (double)SYNTHETIC_LATITUDE_CELLS;
    double latitude_upper = -0.5 * pi
        + pi * (double)(latitude_cell + 1)
        / (double)SYNTHETIC_LATITUDE_CELLS;

    /* Exact unit-sphere area of one longitude-latitude cell. */
    return (sin(latitude_upper) - sin(latitude_lower))
        * (2.0 * pi / (double)SYNTHETIC_LONGITUDE_CELLS);
}

static double manufactured_perturbation(int global_cell, int radial_node,
                                        int radial_nodes)
{
    int latitude_cell = global_cell / SYNTHETIC_LONGITUDE_CELLS;
    int longitude_cell = global_cell % SYNTHETIC_LONGITUDE_CELLS;
    double latitude = -0.5 * pi
        + pi * ((double)latitude_cell + 0.5)
        / (double)SYNTHETIC_LATITUDE_CELLS;
    double longitude = 2.0 * pi
        * ((double)longitude_cell + 0.5)
        / (double)SYNTHETIC_LONGITUDE_CELLS;
    double horizontal = sin(2.0 * longitude) * cos(latitude)
        + 0.35 * cos(3.0 * longitude) * sin(2.0 * latitude);
    double radial = cos(3.0 * pi * (double)radial_node
                        / (double)(radial_nodes - 1));

    return 0.03 * horizontal * radial;
}

static void fill_fields(const double *rho, const double *gravity,
                        const double *tref, const double *alpha,
                        int global_start, int local_cells, int radial_nodes,
                        double atemp, int production_scaled,
                        double *absolute, double *anomaly)
{
    int cell, node;

    for (cell = 0; cell < local_cells; ++cell) {
        int global_cell = global_start + cell;
        for (node = 0; node < radial_nodes; ++node) {
            double perturbation = manufactured_perturbation(
                global_cell, node, radial_nodes);
            if (production_scaled) {
                anomaly[cell * radial_nodes + node] =
                    atemp * rho[node] * alpha[node] * perturbation
                    * gravity[node];
                absolute[cell * radial_nodes + node] =
                    atemp * rho[node] * alpha[node]
                    * (tref[node] + perturbation) * gravity[node];
            }
            else {
                anomaly[cell * radial_nodes + node] =
                    rho[node] * alpha[node] * perturbation;
                absolute[cell * radial_nodes + node] =
                    rho[node] * alpha[node] * (tref[node] + perturbation);
            }
        }
    }
}

static void fill_reference_field(const double *rho, const double *gravity,
                                 const double *tref, const double *alpha,
                                 int local_cells, int radial_nodes,
                                 double atemp, int production_scaled,
                                 double *field)
{
    int cell, node;

    for (cell = 0; cell < local_cells; ++cell)
        for (node = 0; node < radial_nodes; ++node) {
            if (production_scaled)
                field[cell * radial_nodes + node] =
                    atemp * rho[node] * alpha[node] * tref[node]
                    * gravity[node];
            else
                field[cell * radial_nodes + node] =
                    rho[node] * alpha[node] * tref[node];
        }
}

static void project_weighted(const double *field, const double *weights,
                             int local_cells, int radial_nodes,
                             MPI_Comm communicator, double *projected)
{
    double *local_numerator;
    double *global_numerator;
    double local_weight = 0.0;
    double global_weight;
    int cell, node;

    local_numerator = (double *)calloc((size_t)radial_nodes, sizeof(double));
    global_numerator = (double *)calloc((size_t)radial_nodes, sizeof(double));
    if (local_numerator == NULL || global_numerator == NULL)
        MPI_Abort(communicator, 2);

    for (cell = 0; cell < local_cells; ++cell) {
        local_weight += weights[cell];
        for (node = 0; node < radial_nodes; ++node)
            local_numerator[node] +=
                weights[cell] * field[cell * radial_nodes + node];
    }

    MPI_Allreduce(local_numerator, global_numerator, radial_nodes, MPI_DOUBLE,
                  MPI_SUM, communicator);
    MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM,
                  communicator);
    if (!isfinite(global_weight) || global_weight <= 0.0)
        MPI_Abort(communicator, 2);

    for (cell = 0; cell < local_cells; ++cell)
        for (node = 0; node < radial_nodes; ++node)
            projected[cell * radial_nodes + node] =
                field[cell * radial_nodes + node]
                - global_numerator[node] / global_weight;

    free(global_numerator);
    free(local_numerator);
}

static struct norms global_weighted_norms(const double *field,
                                          const double *weights,
                                          int local_cells, int radial_nodes,
                                          MPI_Comm communicator)
{
    struct norms result;
    double local_square = 0.0;
    double global_square;
    double local_maximum = 0.0;
    double global_maximum;
    double local_weight = 0.0;
    double global_weight;
    int local_bad = 0;
    int global_bad;
    int cell, node;

    for (cell = 0; cell < local_cells; ++cell) {
        if (!isfinite(weights[cell]) || weights[cell] <= 0.0)
            local_bad = 1;
        for (node = 0; node < radial_nodes; ++node) {
            double value = field[cell * radial_nodes + node];
            if (!isfinite(value)) {
                local_bad = 1;
                continue;
            }
            local_square += weights[cell] * value * value;
            if (fabs(value) > local_maximum)
                local_maximum = fabs(value);
            local_weight += weights[cell];
        }
    }

    MPI_Allreduce(&local_square, &global_square, 1, MPI_DOUBLE, MPI_SUM,
                  communicator);
    MPI_Allreduce(&local_maximum, &global_maximum, 1, MPI_DOUBLE, MPI_MAX,
                  communicator);
    MPI_Allreduce(&local_weight, &global_weight, 1, MPI_DOUBLE, MPI_SUM,
                  communicator);
    MPI_Allreduce(&local_bad, &global_bad, 1, MPI_INT, MPI_MAX, communicator);

    result.weighted_l2 = sqrt(global_square);
    result.weighted_rms = global_weight > 0.0
        ? sqrt(global_square / global_weight) : NAN;
    result.linf = global_maximum;
    result.total_weight = global_weight;
    result.finite = !global_bad && isfinite(result.weighted_l2)
        && isfinite(result.weighted_rms) && isfinite(result.linf)
        && isfinite(result.total_weight) && result.total_weight > 0.0;
    return result;
}

static void serial_weighted_means(const double *rho, const double *tref,
                                  const double *alpha, int radial_nodes,
                                  double *absolute_mean,
                                  double *anomaly_mean)
{
    double total_weight = 0.0;
    int global_cell, node;

    for (node = 0; node < radial_nodes; ++node) {
        absolute_mean[node] = 0.0;
        anomaly_mean[node] = 0.0;
    }
    for (global_cell = 0; global_cell < SYNTHETIC_SURFACE_CELLS;
         ++global_cell) {
        double weight = surface_weight(global_cell);
        total_weight += weight;
        for (node = 0; node < radial_nodes; ++node) {
            double perturbation = manufactured_perturbation(
                global_cell, node, radial_nodes);
            double coefficient = rho[node] * alpha[node];
            anomaly_mean[node] += weight * coefficient * perturbation;
            absolute_mean[node] +=
                weight * coefficient * (tref[node] + perturbation);
        }
    }
    for (node = 0; node < radial_nodes; ++node) {
        absolute_mean[node] /= total_weight;
        anomaly_mean[node] /= total_weight;
    }
}

static int exceeds_absolute_tolerance(struct norms value, double tolerance)
{
    return !value.finite || value.weighted_l2 > tolerance
        || value.weighted_rms > tolerance || value.linf > tolerance;
}

static double maximum_relative_norm(struct norms residual,
                                    struct norms reference)
{
    double l2, rms, linf, maximum;

    if (!residual.finite || !reference.finite
        || reference.weighted_l2 <= DBL_MIN
        || reference.weighted_rms <= DBL_MIN
        || reference.linf <= DBL_MIN)
        return INFINITY;
    l2 = residual.weighted_l2 / reference.weighted_l2;
    rms = residual.weighted_rms / reference.weighted_rms;
    linf = residual.linf / reference.linf;
    maximum = l2 > rms ? l2 : rms;
    return maximum > linf ? maximum : linf;
}

int main(int argc, char **argv)
{
    double rho[MAX_RADIAL_NODES];
    double gravity[MAX_RADIAL_NODES];
    double tref[MAX_RADIAL_NODES];
    double alpha[MAX_RADIAL_NODES];
    double requested_atemp = 1.0;
    double effective_atemp = 1.0;
    const char *atemp_source = "UNIT_DEFAULT";
    double *weights;
    double *absolute;
    double *anomaly;
    double *projected_absolute;
    double *projected_anomaly;
    double *difference;
    double *serial_absolute_mean;
    double *serial_anomaly_mean;
    struct norms coefficient_identity;
    struct norms coefficient_zero;
    struct norms coefficient_signal;
    struct norms decomposition_absolute;
    struct norms decomposition_anomaly;
    struct norms scaled_identity;
    struct norms scaled_zero;
    struct norms scaled_reference;
    double scaled_identity_relative;
    double scaled_zero_relative;
    int rank, ranks, radial, radial_min, radial_max;
    int global_start, local_cells, allocation_count;
    int cell, node;
    int failed;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);

    if (argc < 2 || argc > 3
        || (argc == 3
            && !parse_atemp(argv[2], &requested_atemp, &effective_atemp))) {
        if (rank == 0)
            fprintf(stderr, "usage: %s refstate_ALA_strict.txt [Atemp]\n",
                    argv[0]);
        MPI_Finalize();
        return 2;
    }
    if (argc == 3)
        atemp_source = "CLI";

    radial = load_reference(argv[1], rho, gravity, tref, alpha);
    MPI_Allreduce(&radial, &radial_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&radial, &radial_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (radial_min < 2 || radial_min != radial_max) {
        if (rank == 0)
            fprintf(stderr,
                    "failed to load identical strict refstate on all ranks: "
                    "minimum_code=%d maximum_code=%d\n",
                    radial_min, radial_max);
        MPI_Finalize();
        return 2;
    }
    radial = radial_min;

    partition_surface(rank, ranks, &global_start, &local_cells);
    allocation_count = local_cells > 0 ? local_cells * radial : 1;
    weights = (double *)malloc((size_t)(local_cells > 0 ? local_cells : 1)
                               * sizeof(double));
    absolute = (double *)malloc((size_t)allocation_count * sizeof(double));
    anomaly = (double *)malloc((size_t)allocation_count * sizeof(double));
    projected_absolute = (double *)malloc(
        (size_t)allocation_count * sizeof(double));
    projected_anomaly = (double *)malloc(
        (size_t)allocation_count * sizeof(double));
    difference = (double *)malloc((size_t)allocation_count * sizeof(double));
    serial_absolute_mean = (double *)malloc((size_t)radial * sizeof(double));
    serial_anomaly_mean = (double *)malloc((size_t)radial * sizeof(double));
    if (weights == NULL || absolute == NULL || anomaly == NULL
        || projected_absolute == NULL || projected_anomaly == NULL
        || difference == NULL || serial_absolute_mean == NULL
        || serial_anomaly_mean == NULL)
        MPI_Abort(MPI_COMM_WORLD, 2);

    for (cell = 0; cell < local_cells; ++cell)
        weights[cell] = surface_weight(global_start + cell);

    fill_fields(rho, gravity, tref, alpha, global_start, local_cells, radial,
                effective_atemp, 0, absolute, anomaly);
    project_weighted(absolute, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_absolute);
    project_weighted(anomaly, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_anomaly);
    for (node = 0; node < local_cells * radial; ++node)
        difference[node] = projected_absolute[node]
            - projected_anomaly[node];
    coefficient_identity = global_weighted_norms(
        difference, weights, local_cells, radial, MPI_COMM_WORLD);
    coefficient_signal = global_weighted_norms(
        projected_anomaly, weights, local_cells, radial, MPI_COMM_WORLD);

    /* Compare the distributed result with a fixed-order serial oracle. */
    serial_weighted_means(rho, tref, alpha, radial,
                          serial_absolute_mean, serial_anomaly_mean);
    for (cell = 0; cell < local_cells; ++cell)
        for (node = 0; node < radial; ++node)
            difference[cell * radial + node] =
                projected_absolute[cell * radial + node]
                - (absolute[cell * radial + node]
                   - serial_absolute_mean[node]);
    decomposition_absolute = global_weighted_norms(
        difference, weights, local_cells, radial, MPI_COMM_WORLD);
    for (cell = 0; cell < local_cells; ++cell)
        for (node = 0; node < radial; ++node)
            difference[cell * radial + node] =
                projected_anomaly[cell * radial + node]
                - (anomaly[cell * radial + node]
                   - serial_anomaly_mean[node]);
    decomposition_anomaly = global_weighted_norms(
        difference, weights, local_cells, radial, MPI_COMM_WORLD);

    fill_reference_field(rho, gravity, tref, alpha, local_cells, radial,
                         effective_atemp, 0, absolute);
    project_weighted(absolute, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_absolute);
    coefficient_zero = global_weighted_norms(
        projected_absolute, weights, local_cells, radial, MPI_COMM_WORLD);

    /*
     * Diagnostic only: mimic production thermal force scaling and operation
     * order, Atemp*rho*alpha*T*g, while retaining the synthetic projector.
     * Raw absolute cancellation is not gated at 1e-10; normalized residuals
     * are, because the unprojected force is O(1e8) for this configuration.
     */
    fill_fields(rho, gravity, tref, alpha, global_start, local_cells, radial,
                effective_atemp, 1, absolute, anomaly);
    project_weighted(absolute, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_absolute);
    project_weighted(anomaly, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_anomaly);
    for (node = 0; node < local_cells * radial; ++node)
        difference[node] = projected_absolute[node]
            - projected_anomaly[node];
    scaled_identity = global_weighted_norms(
        difference, weights, local_cells, radial, MPI_COMM_WORLD);

    fill_reference_field(rho, gravity, tref, alpha, local_cells, radial,
                         effective_atemp, 1, absolute);
    scaled_reference = global_weighted_norms(
        absolute, weights, local_cells, radial, MPI_COMM_WORLD);
    project_weighted(absolute, weights, local_cells, radial, MPI_COMM_WORLD,
                     projected_absolute);
    scaled_zero = global_weighted_norms(
        projected_absolute, weights, local_cells, radial, MPI_COMM_WORLD);
    scaled_identity_relative = maximum_relative_norm(
        scaled_identity, scaled_reference);
    scaled_zero_relative = maximum_relative_norm(scaled_zero, scaled_reference);

    failed = exceeds_absolute_tolerance(coefficient_identity,
                                        coefficient_tolerance)
        || exceeds_absolute_tolerance(coefficient_zero,
                                      coefficient_tolerance)
        || exceeds_absolute_tolerance(decomposition_absolute,
                                      coefficient_tolerance)
        || exceeds_absolute_tolerance(decomposition_anomaly,
                                      coefficient_tolerance)
        || !coefficient_signal.finite
        || coefficient_signal.linf <= minimum_signal
        || !isfinite(scaled_identity_relative)
        || !isfinite(scaled_zero_relative)
        || scaled_identity_relative > scaled_relative_tolerance
        || scaled_zero_relative > scaled_relative_tolerance;

    if (rank == 0) {
        printf("MPI_THERMODYNAMIC_CLOSURE "
               "scope=STANDALONE_SYNTHETIC_THERMAL_ONLY "
               "production_projection=NOT_TESTED phase=NOT_TESTED "
               "composition=NOT_TESTED ranks=%d global_surface_cells=%d "
               "radial=%d projector=SYNTHETIC_SPHERICAL_CELL_AREA\n",
               ranks, SYNTHETIC_SURFACE_CELLS, radial);
        printf("COEFFICIENT_IDENTITY weighted_L2=%.17e weighted_RMS=%.17e "
               "Linf=%.17e tolerance=%.1e\n",
               coefficient_identity.weighted_l2,
               coefficient_identity.weighted_rms,
               coefficient_identity.linf, coefficient_tolerance);
        printf("COEFFICIENT_ZERO_REFERENCE weighted_L2=%.17e "
               "weighted_RMS=%.17e Linf=%.17e tolerance=%.1e\n",
               coefficient_zero.weighted_l2,
               coefficient_zero.weighted_rms,
               coefficient_zero.linf, coefficient_tolerance);
        printf("DECOMPOSITION_ORACLE absolute_weighted_L2=%.17e "
               "absolute_weighted_RMS=%.17e absolute_Linf=%.17e "
               "anomaly_weighted_L2=%.17e anomaly_weighted_RMS=%.17e "
               "anomaly_Linf=%.17e signal_Linf=%.17e tolerance=%.1e\n",
               decomposition_absolute.weighted_l2,
               decomposition_absolute.weighted_rms,
               decomposition_absolute.linf,
               decomposition_anomaly.weighted_l2,
               decomposition_anomaly.weighted_rms,
               decomposition_anomaly.linf, coefficient_signal.linf,
               coefficient_tolerance);
        printf("PRODUCTION_SCALED_DIAGNOSTIC Atemp_source=%s "
               "gravity_source=REFSTATE_COLUMN_2 requested_Atemp=%.17e "
               "effective_float_Atemp=%.17e "
               "identity_raw_weighted_L2=%.17e "
               "identity_raw_weighted_RMS=%.17e identity_raw_Linf=%.17e "
               "identity_max_relative=%.17e zero_raw_weighted_L2=%.17e "
               "zero_raw_weighted_RMS=%.17e zero_raw_Linf=%.17e "
               "zero_max_relative=%.17e relative_tolerance=%.1e\n",
               atemp_source, requested_atemp, effective_atemp,
               scaled_identity.weighted_l2, scaled_identity.weighted_rms,
               scaled_identity.linf, scaled_identity_relative,
               scaled_zero.weighted_l2, scaled_zero.weighted_rms,
               scaled_zero.linf, scaled_zero_relative,
               scaled_relative_tolerance);
        printf("MPI_THERMODYNAMIC_CLOSURE status=%s\n",
               failed ? "FAIL" : "PASS");
        fflush(stdout);
    }

    free(serial_anomaly_mean);
    free(serial_absolute_mean);
    free(difference);
    free(projected_anomaly);
    free(projected_absolute);
    free(anomaly);
    free(absolute);
    free(weights);
    MPI_Finalize();
    return failed ? 1 : 0;
}
