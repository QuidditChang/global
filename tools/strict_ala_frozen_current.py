#!/usr/bin/env python3
"""Prepare and assess a fresh, bounded strict-ALA rheol7 frozen experiment.

No historical POD, hash or convergence trajectory is treated as current data.
Only configuration/provenance/analysis: this tool never modifies solver physics.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

CASES = ("ADJOINT", "BASE", "TIGHT", "BPI", "UNSCALED")
VS = "CitcomS.solver.vsolver"
COMMON = {
    ("CitcomS", "steps"): "1",
    (VS, "piterations"): "60",
    (VS, "ala_pcg_restart_interval"): "50",
    (VS, "ala_outer_solver"): "fgmres",
    (VS, "ala_stage_abc_production_logging"): "on",
    (VS, "ala_stage_abc_adjoint_diagnostic"): "off",
    (VS, "ala_stage_e_diagnostic"): "on",
    (VS, "ala_schur_diagnostic"): "off",
    (VS, "ala_schur_diagnostic_only"): "off",
    (VS, "ala_schur_diagnostic_resume"): "off",
    (VS, "ala_depth_diagnostics"): "on",
    (VS, "ala_depth_diagnostic_interval"): "5",
    (VS, "ala_coarse_residual_diagnostics"): "on",
    (VS, "ala_coarse_residual_interval"): "5",
    (VS, "ala_viscosity_spectrum_diagnostics"): "on",
    ("CitcomS.controller", "profileMonitoringFrequency"): "1",
}
VARIANTS = {
    "ADJOINT": {(VS, "ala_stage_abc_adjoint_diagnostic"): "on",
                (VS, "ala_stage_e_diagnostic"): "off",
                (VS, "ala_schur_diagnostic_tight_tolerance"): "1e-8"},
    "BASE": {},
    "TIGHT": {(VS, "ala_inner_accuracy_max"): "1e-4",
              (VS, "ala_inner_accuracy_factor"): "1e-4"},
    "BPI": {(VS, "ala_shallow_patch_preconditioner"): "off"},
    "UNSCALED": {(VS, "ala_shallow_patch_mid_action_scale"): "1.0",
                 (VS, "ala_shallow_patch_transition_action_scale"): "1.0"},
}
INPUTS = ("cmbhf_ALA_strict.cfg", "GLB.coor.global.dat",
          "refstate_ALA_strict.txt", "interval_ALA_strict.txt")


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                    allow_nan=False) + "\n")


def cfg_read(path):
    result, section = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
        elif "=" in line:
            key, value = (s.strip() for s in line.split("=", 1))
            if (section, key) in result:
                raise ValueError("duplicate cfg key: " + str((section, key)))
            result[section, key] = value
    return result


def cfg_write(path, values):
    sections = dict.fromkeys(s for s, _ in values)
    text = "# Generated frozen experiment; canonical configuration is unchanged.\n"
    for section in sections:
        text += "\n[" + section + "]\n"
        text += "".join(k + " = " + v + "\n" for (s, k), v in values.items()
                        if s == section)
    Path(path).write_text(text)


def binary_files(code):
    files = {code / "bin/pycitcoms"}
    for pattern in ("module/**/*.so", "python/**/*.so", "lib/**/*.so", "lib/**/libCitcomS.a",
                    "module/**/libCitcomSLibmodule.a"):
        files.update(p for p in code.glob(pattern) if p.is_file())
    if not all(p.is_file() for p in files):
        raise ValueError("missing built bin/pycitcoms")
    return {str(p.resolve()): digest(p) for p in sorted(files)}


def git_head(path):
    return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"],
                                   text=True).strip()


def stamp_build(code):
    # Run this only after config_script completed successfully, as in the README.
    subprocess.run(["git", "-C", str(code), "diff", "--exit-code", "HEAD",
                    "--", "lib", "CitcomS", "module"], check=True)
    write_json(code / "frozen_current_build.json",
               {"code_commit": git_head(code), "binaries": binary_files(code)})


def validate_inputs(runs):
    import numpy as np
    cfg = cfg_read(runs / INPUTS[0])
    expected = {(VS, "compressible_formulation"): "ala",
                (VS, "uzawa"): "ala_cg", (VS, "remove_rigid_rotation"): "0",
                (VS, "ala_beta_element_source"): "interval",
                (VS, "ala_inner_accuracy_max"): "1e-2",
                (VS, "ala_shallow_patch_preconditioner"): "on",
                ("CitcomS.solver.visc", "rheol"): "7",
                ("CitcomS.solver.visc", "SDEPV"): "off",
                ("CitcomS.solver.visc", "PDEPV"): "off",
                ("CitcomS.solver.param", "refstate_file"): INPUTS[2],
                ("CitcomS.solver.param", "ala_beta_interval_file"): INPUTS[3]}
    for key, value in expected.items():
        if cfg.get(key) != value:
            raise ValueError("review changed launch contract: " + str(key))
    mesh = np.loadtxt(runs / INPUTS[1], skiprows=1)
    ref = np.loadtxt(runs / INPUTS[2])
    interval = np.loadtxt(runs / INPUTS[3])
    if mesh.shape != (65, 2) or ref.shape != (65, 8) or interval.shape != (64, 4):
        raise ValueError("unexpected mesh/reference/interval dimensions")
    if not all(np.isfinite(x).all() for x in (mesh, ref, interval)):
        raise ValueError("nonfinite input")
    if not (ref > 0).all() or not (interval[:, 3] > 0).all():
        raise ValueError("nonpositive reference coefficient")
    if not np.allclose(interval[:, 1:3], np.column_stack((mesh[:-1, 1], mesh[1:, 1])),
                       rtol=0, atol=1e-12):
        raise ValueError("interval/mesh mismatch")
    di = float(cfg["CitcomS.solver", "dissipation_number"])
    if not np.allclose(ref[:, 6], di*ref[:, 3]*ref[:, 1]/ref[:, 4]/ref[:, 5],
                       rtol=1e-10, atol=0):
        raise ValueError("Gamma closure mismatch")
    ms = "CitcomS.solver.mesher"
    if tuple(int(cfg[ms, k]) for k in ("nprocx", "nprocy", "nprocz")) != (4,4,2):
        raise ValueError("experiment requires the current 384-rank layout")
    ic = "CitcomS.solver.ic"
    if cfg.get((ic,"restart"),"off") not in ("off","0"):
        raise ValueError("experiment requires a fresh step-0 state")
    return cfg


def prepare(runs, code, root, hpc):
    cfg = validate_inputs(runs)
    if root.exists():
        raise ValueError("use a fresh experiment directory: " + str(root))
    extra = {}
    if hpc:
        receipt = json.loads((code / "frozen_current_build.json").read_text())
        if receipt != {"code_commit": git_head(code), "binaries": binary_files(code)}:
            raise ValueError("build receipt mismatch: rebuild and stamp")
        # Record every reconstruction file selected around this fixed start age.
        age = float(cfg["CitcomS.solver.param", "start_age"])
        for key in ("vel_bound_file", "lith_age_file", "flag_depth_new_file",
                    "flag_depth_file", "tf_file"):
            prefix = Path(cfg["CitcomS.solver.param", key])
            for year in {math.floor(age), math.ceil(age)}:
                files = sorted(prefix.parent.glob(prefix.name + str(year) + ".*"))
                if not files:
                    raise ValueError("missing reconstruction prefix: " + str(prefix) + str(year))
                extra.update({str(p): digest(p) for p in files if p.is_file()})
    root.mkdir(parents=True)
    manifest = {"schema": 1, "code_commit": git_head(code), "code_directory": str(code),
                "runs_commit": git_head(runs), "cases": {},
                "inputs": {name: digest(runs/name) for name in INPUTS},
                "nuref_header_sha256": digest(code/"lib/Steinberger_nuref.h"),
                "reconstruction_files": extra, "hpc_preflight": hpc}
    if hpc:
        manifest["build"] = receipt
    for case in CASES:
        directory = root/case
        (directory/"DATA/0").mkdir(parents=True)
        (directory/"Restart").mkdir()
        values = {**cfg, **COMMON, **VARIANTS[case]}
        values["CitcomS.solver", "datadir"] = str(directory/"DATA/%RANK")
        values["CitcomS.solver", "datadir_old"] = str(directory/"Restart")
        cfg_write(directory/"case.cfg", values)
        for name in INPUTS[1:]:
            shutil.copy2(runs/name, directory/name)
        manifest["cases"][case] = {
            "cfg_sha256": digest(directory/"case.cfg"),
            "changes": {s+"."+k: [cfg.get((s,k)),v]
                        for (s,k),v in values.items() if cfg.get((s,k)) != v}}
    write_json(root/"manifest.json", manifest)
    print(root)


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def analyze(root):
    manifest = json.loads((root/"manifest.json").read_text())
    result = {"valid": True, "errors": [], "cases": {}, "interpretation": []}
    if manifest.get("hpc_preflight"):
        try:
            if binary_files(Path(manifest["code_directory"])) != manifest["build"]["binaries"]:
                raise ValueError("binaries changed during campaign")
            if any(digest(p)!=h for p,h in manifest["reconstruction_files"].items()):
                raise ValueError("reconstruction inputs changed during campaign")
        except (ValueError,OSError) as exc:
            result["valid"] = False
            result["errors"].append(str(exc))
    hashes = []
    for case in CASES:
        d = root/case
        try:
            if digest(d/"case.cfg") != manifest["cases"][case]["cfg_sha256"]:
                raise ValueError("cfg changed after preflight")
            for name in INPUTS[1:]:
                if digest(d/name) != manifest["inputs"][name]:
                    raise ValueError("input changed after preflight: " + name)
            if int((d/"exit_status.txt").read_text()) != 0:
                raise ValueError("nonzero application exit status")
            logs = list((d/"DATA/0").glob("global_AhatP*.log"))
            if len(logs) != 1:
                raise ValueError("expected exactly one solver log")
            text = logs[0].read_text(errors="replace")
            frozen = re.findall(r"STRICT_ALA_FROZEN_CURRENT phase=(begin|end) physical=([0-9a-f]+) viscosity=([0-9a-f]+) matched=1 step=0 time=([^\s]+)", text)
            if (len(frozen)!=2 or frozen[0][0]!="begin" or frozen[1][0]!="end"
                    or frozen[0][1:]!=frozen[1][1:]):
                raise ValueError("missing/failed physical and multilevel-viscosity freeze guard")
            hashes.append(tuple(frozen[0][1:]))
            if case == "ADJOINT":
                adj = rows(d/"DATA/0/global.strict_ala_stage_B_adjoint.csv")
                gauge = rows(d/"DATA/0/global.strict_ala_stage_B_gauge.csv")
                if len(adj) != 42 or len(gauge) != 4:
                    # Seven pressure probes x two velocities x D/C/B.
                    raise ValueError("incomplete operator audit rows")
                for records,keys in ((adj,("delta_scale","repeatability_relative")),
                                     (gauge,("achieved_first","achieved_second",
                                             "qT_S_gamma_q","BT_split_relative",
                                             "repeat_action_relative"))):
                    if not all(math.isfinite(float(r[k])) for r in records for k in keys):
                        raise ValueError("nonfinite operator audit sample")
                defect = max(float(r["delta_scale"]) for r in adj)
                repeat = max(float(r["repeatability_relative"]) for r in adj)
                achieved = max(float(r[k]) for r in gauge
                               for k in ("achieved_first", "achieved_second"))
                curvature = min(float(r["qT_S_gamma_q"]) for r in gauge)
                split = max(float(r["BT_split_relative"]) for r in gauge)
                action_repeat = max(float(r["repeat_action_relative"]) for r in gauge)
                result["cases"][case] = dict(adjoint_defect=defect, repeatability=repeat,
                    achieved_inner=achieved, minimum_sampled_curvature=curvature,
                    transpose_split_defect=split, action_repeatability=action_repeat,
                    operator_gate_pass=(defect<=1e-8 and repeat<=1e-8 and split<=1e-8
                                        and action_repeat<=1e-6 and achieved<=1.01e-8 and curvature>0))
                if not all(math.isfinite(x) for x in (defect,repeat,achieved,curvature)):
                    raise ValueError("nonfinite operator audit")
                if defect > 1e-8 or repeat > 1e-8:
                    result["interpretation"].append("OPERATOR_OR_MPI_ADJOINT_FAILURE")
                if achieved > 1.01e-8:
                    result["interpretation"].append("TIGHT_VELOCITY_INVERSE_NOT_VERIFIED")
                if curvature <= 0:
                    result["interpretation"].append("NONPOSITIVE_SAMPLED_SCHUR_CURVATURE")
                if split>1e-8 or action_repeat>1e-6:
                    result["interpretation"].append("TRANSPOSE_SPLIT_OR_ACTION_REPEATABILITY_FAILURE")
                continue
            it = rows(d/"DATA/0/global.strict_ala_stage_C_iterations.csv")
            inner = rows(d/"DATA/0/global.strict_ala_stage_C_inner_solves.csv")
            if not it or not inner or it[-1]["final_iterate"] != "1":
                raise ValueError("incomplete trajectory")
            for r in it:
                if not all(math.isfinite(float(r[k])) for k in
                           ("continuity_relative","momentum_relative","krylov_drift")):
                    raise ValueError("nonfinite trajectory")
            cfg = cfg_read(d/"case.cfg")
            accepted = (float(it[-1]["continuity_relative"]) <= float(cfg[VS,"tole_compressibility"])
                        and float(it[-1]["momentum_relative"]) <= float(cfg[VS,"ala_unaugmented_momentum_tolerance"]))
            if not accepted and int(it[-1]["iteration"]) != 60:
                raise ValueError("stopped early without joint acceptance")
            summary = {"joint_converged": bool(accepted),
                "final_iteration": int(it[-1]["iteration"]),
                "best_continuity": min(float(r["continuity_relative"]) for r in it),
                "final_continuity": float(it[-1]["continuity_relative"]),
                "final_momentum": float(it[-1]["momentum_relative"]),
                "inner_failures": sum(r["status"]!="CONVERGED" for r in inner),
                "inner_seconds": sum(float(r["seconds"]) for r in inner),
                "trajectory": {r["iteration"]: float(r["continuity_relative"]) for r in it},
                "physical_hash": frozen[0][1], "viscosity_hash": frozen[0][2]}
            result["cases"][case] = summary
            if summary["inner_failures"]:
                result["interpretation"].append(case+":INNER_SOLVE_TARGET_MISSED")
        except (ValueError, OSError, KeyError, IndexError) as exc:
            result["valid"] = False
            result["errors"].append(case+": "+str(exc))
    if len(set(hashes)) > 1:
        result["valid"] = False
        result["errors"].append("different initialized physics/viscosity across cases")
    if result["valid"]:
        base = result["cases"]["BASE"]
        for case in ("TIGHT","BPI","UNSCALED"):
            other = result["cases"][case]
            common = set(base["trajectory"]) & set(other["trajectory"])
            last = max(common, key=int)
            result["cases"][case]["comparison_iteration"] = int(last)
            result["cases"][case]["continuity_ratio_to_base"] = (
                other["trajectory"][last]/max(base["trajectory"][last],1e-300))
        result["interpretation"].append("COMPARE_MATCHED_ITERATIONS_AND_COST; ratios are evidence, not causal proof")
        result["interpretation"].append("BASE_CONVERGED" if base["joint_converged"] else
            "BASE_NOT_CONVERGED_WITHIN_60; inspect fresh depth/restart/Hessenberg diagnostics")
    result["solver_attribution_allowed"] = (result["valid"]
        and result["cases"].get("ADJOINT",{}).get("operator_gate_pass",False)
        and all(result["cases"].get(c,{}).get("inner_failures",1)==0
                for c in ("BASE","TIGHT","BPI","UNSCALED")))
    write_json(root/"analysis.json", result)
    print(json.dumps(result, indent=2))
    return 0 if result["valid"] else 2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="action", required=True)
    pre = sub.add_parser("prepare")
    pre.add_argument("--runs", type=Path, required=True)
    pre.add_argument("--code", type=Path, required=True)
    pre.add_argument("--root", type=Path, required=True)
    pre.add_argument("--hpc", action="store_true")
    stamp = sub.add_parser("stamp-build")
    stamp.add_argument("--code", type=Path, required=True)
    ana = sub.add_parser("analyze")
    ana.add_argument("--root", type=Path, required=True)
    a = p.parse_args()
    if a.action == "prepare":
        prepare(a.runs.resolve(), a.code.resolve(), a.root.resolve(), a.hpc)
    elif a.action == "stamp-build":
        stamp_build(a.code.resolve())
    else:
        return analyze(a.root.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
