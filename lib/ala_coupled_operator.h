#ifndef ALA_COUPLED_OPERATOR_H
#define ALA_COUPLED_OPERATOR_H

struct All_variables;

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

#endif
