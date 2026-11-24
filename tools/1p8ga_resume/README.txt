1.8 Ga resume: original producer 0608be89a18479a0d22374eb9b395be940c49d86

Only the Full_tracer_advection.c hunks from b7ace82 and 4c2c781 are
backported. Pices_checkpoint.c permits the explicitly named original
producer while requiring exact equality of every other metadata byte.
No rheology, thermal, Stokes, CFL or checkpoint layout changes.

Submit submit.lsf with CODE_DIR, SCRIPTDIR and RUNROOT in bsub -env.
SCRIPTDIR must be a clean checkout at 8ef9049d3fb30147d8c4fd07d21f20cd4272b5f6.
CODE_DIR must be the built codex/pices-1p8ga-host-fix checkout.
RUNROOT is the original stable run_links symlink, not a dated scratch path.

Preparation validates original Runtime, all 384 checkpoints at step 39400,
collective manifests and payload/state hashes. It copies the selected
checkpoint to Restart/host_fix_39400 before altering input.cfg. It refuses
an existing backup (inspect an interrupted attempt; do not blindly delete).
The new Runtime_host_fix_<commit> is frozen separately. Original Runtime
and solver_commit.txt remain provenance for the original run; the new
restart_manifest.json and new Runtime manifest record the resumed solver.

400 slots: 384 solver plus 16 postprocessing; seven-day rollover retained.
300000 step ceiling retained. Existing geological-age termination stops
at modern age within the original floating-point time resolution.
No separate 20-step diagnostic job: direct continuation requested by user.
Local build and regression tests do not constitute a replay of the actual
HPC step-39409 failure; production checkpoint preflight remains mandatory.
