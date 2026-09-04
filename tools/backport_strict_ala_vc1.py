#!/usr/bin/env python3
"""Apply the diagnostic-only VC1 hooks to the pinned 2628980 source tree."""

import argparse
import hashlib
import json
import pathlib
import subprocess


DRIVE_SHA256 = "9bb62ec88d6725f9c7b250e3780fa2bfa9650e5803ec66c2b7f402a9530aac9b"
STOKES_SHA256 = "e0d5b0ab89a932449edfdc83a05fe50f45dfad855a29ab95670666b3a5a92f1d"


def digest(data):
    return hashlib.sha256(data.encode()).hexdigest()


def replace_once(text, old, new, label):
    count = text.count(old)
    if count != 1:
        raise SystemExit("VC1 backport context %s expected once, found %d" %
                         (label, count))
    return text.replace(old, new, 1)


def git_text(repo, commit, path):
    return subprocess.check_output(
        ["git", "-C", str(repo), "show", "%s:%s" % (commit, path)],
        text=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-repo", required=True)
    parser.add_argument("--source-tree", required=True)
    parser.add_argument("--diagnostic-commit", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()

    repo = pathlib.Path(args.source_repo)
    tree = pathlib.Path(args.source_tree)
    drive_path = tree / "lib" / "Drive_solvers.c"
    stokes_path = tree / "lib" / "Stokes_flow_Incomp.c"
    drive = drive_path.read_text()
    stokes = stokes_path.read_text()
    if digest(drive) != DRIVE_SHA256 or digest(stokes) != STOKES_SHA256:
        raise SystemExit("pinned 2628980 solver source hash mismatch")

    diagnostic = git_text(repo, args.diagnostic_commit, "lib/Drive_solvers.c")
    helper_start = diagnostic.index("/* VC1 is deliberately process-level continuation:")
    helper_end = diagnostic.index("\n\n/************************************************************/", helper_start)
    helper = diagnostic[helper_start:helper_end]

    drive = replace_once(
        drive, "#include <float.h>\n#include <string.h>",
        "#include <float.h>\n#include <errno.h>\n#include <stdint.h>\n"
        "#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>", "headers")
    drive = replace_once(
        drive, "double global_vdot();\ndouble vnorm_nonnewt();",
        "double global_vdot();\ndouble global_pdot();\ndouble vnorm_nonnewt();\n"
        "void myerror(struct All_variables *,char *);", "declarations")
    drive = replace_once(
        drive, "static void write_stokes_diagnostics(struct All_variables *E);\n",
        "static void write_stokes_diagnostics(struct All_variables *E);\n\n" + helper + "\n",
        "helper block")
    drive = replace_once(
        drive, "  double Udot_mag, dUdot_mag,omega[3];\n  int m,count,i,j,k;",
        "  double Udot_mag, dUdot_mag,omega[3];\n"
        "  double vc1_operator_seconds=0.0,vc1_solver_seconds=0.0,vc1_clock;\n"
        "  unsigned long long vc1_physical_hash_before=0,vc1_physical_hash_after=0;\n"
        "  unsigned long long vc1_local_physical_hash_before=0,vc1_local_physical_hash_after=0;\n"
        "  int vc1_local_physical_mismatch=0,vc1_global_physical_mismatch=0;\n"
        "  int m,count,i,j,k;\n"
        "  const char *vc1_stage=getenv(\"STRICT_ALA_VC1_STAGE\");\n"
        "  const char *vc1_warm_input=getenv(\"STRICT_ALA_VC1_WARM_INPUT\");\n"
        "  const char *vc1_warm_output=getenv(\"STRICT_ALA_VC1_WARM_OUTPUT\");\n"
        "  const char *vc1_compare_left=getenv(\"STRICT_ALA_VC1_COMPARE_LEFT\");\n"
        "  const char *vc1_compare_right=getenv(\"STRICT_ALA_VC1_COMPARE_RIGHT\");",
        "solver locals")
    drive = replace_once(
        drive, "  //velocities_conform_bcs(E,E->U);\n\n  assemble_forces(E,0);",
        "  //velocities_conform_bcs(E,E->U);\n\n"
        "  if(vc1_stage!=NULL && E->monitor.solution_cycles==0 &&\n"
        "     vc1_compare_left!=NULL && vc1_compare_right!=NULL) {\n"
        "    strict_ala_vc1_compare_states(E,vc1_compare_left,vc1_compare_right);\n"
        "    MPI_Barrier(E->parallel.world);\n    MPI_Finalize();\n    exit(EXIT_SUCCESS);\n  }\n"
        "  if(vc1_stage!=NULL && E->monitor.solution_cycles==0 &&\n"
        "     vc1_warm_input!=NULL && vc1_warm_input[0]!='\\0')\n"
        "    strict_ala_vc1_transfer_warm_state(E,vc1_warm_input,0);\n"
        "  if(vc1_stage!=NULL && E->monitor.solution_cycles==0)\n"
        "    vc1_physical_hash_before=strict_ala_vc1_physical_state_hash(\n"
        "        E,&vc1_local_physical_hash_before);\n\n  assemble_forces(E,0);",
        "pre-solve hooks")
    drive = replace_once(
        drive,
        "    get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);\n"
        "    velocities_conform_bcs(E,E->U);\n    construct_stiffness_B_matrix(E);",
        "    vc1_clock=MPI_Wtime();\n"
        "    get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);\n"
        "    velocities_conform_bcs(E,E->U);\n    construct_stiffness_B_matrix(E);\n"
        "    vc1_operator_seconds+=MPI_Wtime()-vc1_clock;", "operator timing")
    drive = replace_once(
        drive, "  }\n\n  solve_constrained_flow_iterative(E);\n\n  if (E->viscosity.SDEPV",
        "  }\n\n"
        "  if(vc1_stage!=NULL && E->monitor.solution_cycles==0) {\n"
        "    strict_ala_vc1_viscosity_audit(E);\n"
        "    if(E->parallel.me==0) {\n"
        "      fprintf(E->fp,\"STRICT_ALA_VC1_FROZEN_STATE_GUARD before=%016llx \"\n"
        "              \"solver_mutable_scope=U_P_only phase_state=derived_from_T_and_cfg\\n\",\n"
        "              vc1_physical_hash_before);\n"
        "      fprintf(E->fp,\"STRICT_ALA_VC1_TIMING operator_rebuild_seconds=%.17e \"\n"
        "              \"fgmres_and_preconditioner_seconds=0.00000000000000000e+00\\n\",\n"
        "              vc1_operator_seconds);\n      fflush(E->fp);\n    }\n  }\n\n"
        "  vc1_clock=MPI_Wtime();\n  solve_constrained_flow_iterative(E);\n"
        "  vc1_solver_seconds+=MPI_Wtime()-vc1_clock;\n\n  if (E->viscosity.SDEPV",
        "first solve")
    drive = replace_once(
        drive,
        "      get_system_viscosity(E,0,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);\n"
        "      velocities_conform_bcs(E,E->U);\n      construct_stiffness_B_matrix(E);\n"
        "      solve_constrained_flow_iterative(E);",
        "      vc1_clock=MPI_Wtime();\n"
        "      get_system_viscosity(E,0,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);\n"
        "      velocities_conform_bcs(E,E->U);\n      construct_stiffness_B_matrix(E);\n"
        "      vc1_operator_seconds+=MPI_Wtime()-vc1_clock;\n"
        "      vc1_clock=MPI_Wtime();\n      solve_constrained_flow_iterative(E);\n"
        "      vc1_solver_seconds+=MPI_Wtime()-vc1_clock;", "nonlinear solve")
    drive = replace_once(
        drive,
        "  if(E->sphere.caps == 12 && E->control.remove_rigid_rotation) {\n"
        "      remove_rigid_rot(E);\n  }\n\n  return;\n}",
        "  if(E->sphere.caps == 12 && E->control.remove_rigid_rotation) {\n"
        "      remove_rigid_rot(E);\n  }\n\n"
        "  if(vc1_stage!=NULL && E->monitor.solution_cycles==0) {\n"
        "    vc1_physical_hash_after=strict_ala_vc1_physical_state_hash(\n"
        "        E,&vc1_local_physical_hash_after);\n"
        "    vc1_local_physical_mismatch=\n"
        "        vc1_local_physical_hash_after!=vc1_local_physical_hash_before;\n"
        "    MPI_Allreduce(&vc1_local_physical_mismatch,&vc1_global_physical_mismatch,\n"
        "                  1,MPI_INT,MPI_MAX,E->parallel.world);\n"
        "    if(vc1_global_physical_mismatch)\n"
        "      myerror(E,\"VC1 frozen T/C/tracer/time state changed during Stokes solve\");\n"
        "    strict_ala_vc1_viscosity_audit(E);\n"
        "    if(E->parallel.me==0) {\n"
        "      fprintf(E->fp,\"STRICT_ALA_VC1_FROZEN_STATE before=%016llx after=%016llx \"\n"
        "              \"bitwise_equal=true phase_state=derived_from_frozen_T_and_cfg\\n\",\n"
        "              vc1_physical_hash_before,vc1_physical_hash_after);\n"
        "      fprintf(E->fp,\"STRICT_ALA_VC1_TIMING operator_rebuild_seconds=%.17e \"\n"
        "              \"fgmres_and_preconditioner_seconds=%.17e\\n\",\n"
        "              vc1_operator_seconds,vc1_solver_seconds);\n      fflush(E->fp);\n    }\n"
        "    if(vc1_warm_output!=NULL && vc1_warm_output[0]!='\\0')\n"
        "      strict_ala_vc1_transfer_warm_state(E,vc1_warm_output,1);\n"
        "  }\n\n  return;\n}", "post-solve hooks")

    stokes = replace_once(
        stokes,
        "                momentum_relative,E->control.ala_unaugmented_momentum_tolerance,\n"
        "                residual_est);\n        fflush(E->fp);",
        "                momentum_relative,E->control.ala_unaugmented_momentum_tolerance,\n"
        "                residual_est);\n"
        "        if(getenv(\"STRICT_ALA_VC1_STAGE\")!=NULL)\n"
        "            fprintf(E->fp,\"STRICT_ALA_VC1_COUPLED_COST \"\n"
        "                    \"K_gamma_rhs_solves=%d velocity_MG_cycles=%d \"\n"
        "                    \"preconditioner_applications=%d source=production_counters\\n\",\n"
        "                    E->control.total_v_solver_calls,\n"
        "                    E->control.total_iteration_cycles,count);\n"
        "        fflush(E->fp);", "cost summary")

    drive_path.write_text(drive)
    stokes_path.write_text(stokes)
    audit = {
        "schema": "strict-ala-vc1-backport-v1",
        "baseline": {"Drive_solvers.c": DRIVE_SHA256,
                     "Stokes_flow_Incomp.c": STOKES_SHA256},
        "diagnostic_commit": args.diagnostic_commit,
        "output": {"Drive_solvers.c": digest(drive),
                   "Stokes_flow_Incomp.c": digest(stokes)},
        "baseline_hashes_verified": True,
        "backport_scope": ["Drive_solvers.c:vc1_state_and_audit_hooks",
                           "Stokes_flow_Incomp.c:read_only_cost_log"],
    }
    pathlib.Path(args.audit).write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
