#ifndef ALA_COUPLED_OPERATOR_H
#define ALA_COUPLED_OPERATOR_H

#include <stddef.h>

struct All_variables;
struct ala_block_vector;

/* Apply the homogeneous strict-ALA saddle-point block
 *
 *        [ K_gamma  (D+C)^T ] [ velocity ]
 *        [  D+C        0    ] [ pressure ] .
 *
 * All fields must be distinct.  The velocity workspace has the same storage
 * contract as the velocity output, including the constrained sentinel at
 * index neq used by the nodal operator.  The pressure output is elementwise.
 */
void apply_ala_coupled_operator(struct All_variables *E,
                                double **velocity, double **pressure,
                                double **momentum, double **continuity,
                                double **velocity_work, int level);

/* Assemble the pressure-independent momentum right-hand side on the finest
 * level.  The legacy force assembly stores f_ext-C^T*p in E->F; this helper
 * removes that historical pressure dependence without changing E->F.
 */
void assemble_ala_pressure_independent_force(struct All_variables *E,
                                             double **force,
                                             double **velocity_work,
                                             int level);

/* Adjacent-level transfer contract for the future coupled V-cycle.  Velocity
 * currently delegates to the established geometric multigrid maps.  P0
 * pressure uses constant prolongation and its exact Euclidean transpose;
 * keeping both here prevents another pressure-only transfer implementation
 * from diverging from the monolithic hierarchy. */
void ala_coupled_prolong_velocity(struct All_variables *E, int coarse_level,
                                  double **coarse, double **fine);
void ala_coupled_restrict_velocity(struct All_variables *E, int fine_level,
                                   double **fine, double **coarse);
void ala_coupled_prolong_pressure_p0(struct All_variables *E,
                                     int coarse_level,
                                     double **coarse, double **fine);
void ala_coupled_restrict_pressure_p0_transpose(struct All_variables *E,
                                                int fine_level,
                                                double **fine,
                                                double **coarse);

/* Stage 9f.0 startup audit.  It checks K_gamma symmetry, the G/G^T adjoint,
 * adjacent velocity and P0-pressure transfer adjoints, pressure mass,
 * duplicate ownership, and the radial coarse-beta restriction on every
 * available level.  It is diagnostic: failures remain visible for the
 * Stage 9f.1/9f.2 implementation rather than silently changing equations. */
void audit_ala_coupled_multilevel_contracts(struct All_variables *E);

/* Stage 9f.1 coupled-patch interface.  A patch owns both fields and solves
 * [K_gamma,P  G_P^T; G_P  -epsilon M_P].  The callback contract lets element
 * and small aggregate kernels share the same V-cycle integration without
 * claiming that the existing velocity-only Vanka is already coupled. */
struct ala_coupled_patch_spec {
    int level;
    int cap;
    const int *elements;
    int element_count;
    double pressure_regularization;
};

struct ala_coupled_patch_stats {
    double pivot_ratio;
    int fallback_count;
    size_t workspace_bytes;
};

typedef int (*ala_coupled_patch_solver)(
    struct All_variables *E, const struct ala_coupled_patch_spec *patch,
    const struct ala_block_vector *residual,
    struct ala_block_vector *correction,
    struct ala_coupled_patch_stats *stats);

#endif
