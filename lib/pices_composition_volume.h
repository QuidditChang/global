#ifndef CITCOMS_PICES_COMPOSITION_VOLUME_H
#define CITCOMS_PICES_COMPOSITION_VOLUME_H
struct All_variables;
#define PICES_VOLUME_FLAVORS 25
struct Pices_volume_summary {
    int initialized, feasible, iterations, hops, step;
    long long surplus_cells, deficit_cells, empty_cells;
    double shell, capacity_min, capacity_max, minimum_transfer_nd;
    double row_error, column_error, species_error, node_error;
    double fraction_error, partition_error, rounding_adjustment_nd;
    double transfer_distance_rms_nd, wall_seconds;
    double target[25], raw[25], corrected[25], node[25], original[25], source[25];
    char reason[160];
};
int pices_composition_volume_enabled(struct All_variables *);
void pices_composition_volume_register(struct All_variables *);
/* initialize does not reconstruct or change existing composition. */
void pices_composition_volume_initialize(struct All_variables *,int restored);
void pices_composition_volume_begin(struct All_variables *);
void pices_composition_volume_transition(struct All_variables *,int cap,int particle,int old_flavor,int new_flavor);
/* Returns feasibility. Shadow failure leaves physics unchanged; active fails. */
int pices_composition_volume_reconstruct(struct All_variables *);
void pices_composition_volume_check(struct All_variables *);
void pices_composition_volume_commit(struct All_variables *);
/* Explicit independent controls only: close files, retain geometry; rebase
 * adopts the already-restored particle weights as the new diagnostic start. */
void pices_composition_volume_reset(struct All_variables *,int rebase);
double pices_composition_volume_weight(struct All_variables *,int cap,int particle);
/* Owned element index is one based, flavor is 0..24. NAN if no feasible map. */
double pices_composition_volume_fraction(struct All_variables *,int element,int flavor);
const struct Pices_volume_summary *pices_composition_volume_latest(struct All_variables *);
#endif
