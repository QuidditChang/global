#ifndef ALA_BLOCK_VECTOR_H
#define ALA_BLOCK_VECTOR_H

struct All_variables;

/* Distributed mixed vector for the strict-ALA saddle-point system.  Velocity
 * uses equation numbering (including the neq sentinel); pressure uses the
 * element/P0 numbering (including the unused zero entry). */
struct ala_block_vector {
    double **velocity;
    double **pressure;
    int level;
};

struct ala_block_vector *ala_block_vector_create(struct All_variables *E,
                                                 int level);
void ala_block_vector_destroy(struct All_variables *E,
                              struct ala_block_vector *vector);
void ala_block_vector_zero(struct All_variables *E,
                           struct ala_block_vector *vector);
void ala_block_vector_copy(struct All_variables *E,
                           const struct ala_block_vector *source,
                           struct ala_block_vector *destination);
void ala_block_vector_scale(struct All_variables *E, double scale,
                            struct ala_block_vector *vector);
void ala_block_vector_axpy(struct All_variables *E, double scale,
                           const struct ala_block_vector *source,
                           struct ala_block_vector *destination);

/* Positive caller-supplied weights define the block Riesz metric.  Velocity
 * uses the duplicate-aware equation dot product and pressure uses its P0
 * volume mass.  Both components are reduced in one MPI collective.  The API
 * intentionally does not remove a pressure mean: with G=D+C, C^T generally
 * makes a constant dynamic pressure observable.  A legacy incompressible
 * pressure-nullspace projection is valid only for an operator that actually
 * has that nullspace (for example the beta=0 limit). */
double ala_block_vector_dot(struct All_variables *E,
                            const struct ala_block_vector *left,
                            const struct ala_block_vector *right,
                            double velocity_weight,
                            double pressure_weight);
double ala_block_vector_norm(struct All_variables *E,
                             const struct ala_block_vector *vector,
                             double velocity_weight,
                             double pressure_weight);

#endif
