#!/usr/bin/env python3
"""Stage A/B/C provenance, cfg, and post-run decision utilities."""
import argparse, csv, hashlib, json, math, os, platform, re, subprocess, sys, time
from pathlib import Path

SCIENCE_KEYS = {
    "ala_shallow_patch_preconditioner",
}
OUTPUT_KEYS = {"datadir", "datafile", "logfile", "datadir_old", "datafile_old"}
REQUIRED_OFF = ("ala_radial_line_preconditioner", "ala_two_level_preconditioner",
                "ala_pressure_multigrid", "ala_pressure_multigrid_galerkin",
                "ala_global_coarse_preconditioner", "ala_geneo_preconditioner")
EXPECTED_PRESSURE_PROBES = ("P1_fixed_random", "P2_radial_smooth",
                            "P5_horizontal_checkerboard", "P6_degree_1",
                            "P6_degree_2", "P14_constant",
                            "P15_density_gauge")
EXPECTED_VELOCITY_PROBES = ("fixed_random", "structured")
EXPECTED_OPERATORS = ("D", "C", "B")
EXPECTED_GAUGE_PROBES = ("P6_degree_1", "P6_degree_2", "P14_constant",
                         "P15_density_gauge")
ADJOINT_NUMERIC_FIELDS = ("lhs", "rhs", "signed_defect", "absolute_defect",
                          "u_norm", "q_norm", "Au_norm", "ATq_norm",
                          "denom_lr", "denom_action", "repeat_lhs_abs",
                          "repeat_rhs_abs", "repeat_ATq_norm", "repeat_Au_norm",
                          "repeatability_relative", "scale_floor", "delta_lr",
                          "delta_scale")
GAUGE_NUMERIC_FIELDS = ("q_norm", "DTq_norm", "CTq_norm", "BTq_norm",
                        "S_gamma_q_norm", "qT_S_gamma_q", "ratio_D", "ratio_B",
                        "BT_sum_norm", "BT_split_absolute", "BT_split_relative",
                        "repeat_action_absolute", "repeat_action_relative",
                        "achieved_first", "achieved_second", "scale_floor")

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def git(repo, *args):
    return subprocess.check_output(["git", "-C", str(repo), *args], text=True,
                                   stderr=subprocess.STDOUT).strip()

def set_key(text, key, value):
    pat = re.compile(r"^(\s*" + re.escape(key) + r"\s*(?:=\s*)?)(\S+)(.*)$", re.M)
    out, n = pat.subn(lambda m: m.group(1) + value + m.group(3), text)
    if not n:
        out += "\n" + key + " = " + value + "\n"
    return out

def cfg_variant(base, out, schwarz):
    text = Path(base).read_text()
    text = set_key(text, "ala_shallow_patch_preconditioner", "on" if schwarz else "off")
    Path(out).write_text(text)

def normalized_cfg_diff(c0, c1):
    def values(path):
        d = {}
        for line in Path(path).read_text().splitlines():
            m = re.match(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.*?)\s*$", line)
            if m and not m.group(1).startswith("#"):
                d[m.group(1)] = m.group(2)
        return d
    a, b = values(c0), values(c1)
    keys = sorted(set(a) | set(b))
    diff = [(k, a.get(k), b.get(k)) for k in keys if a.get(k) != b.get(k)
            and k not in OUTPUT_KEYS]
    return diff

def cfg_values(path):
    values = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r"^\s*([A-Za-z0-9_]+)\s*=\s*(.*?)\s*$", line)
        if m: values[m.group(1)] = m.group(2)
    return values

def cfg_input_files(cfg, missing=None):
    """Resolve cfg keys whose names end in _file, including prefix inputs."""
    cfg = Path(cfg)
    files = []
    for line in cfg.read_text().splitlines():
        m = re.match(r"^\s*([A-Za-z0-9_]+_file)\s*=\s*(\S+)", line)
        if not m or m.group(2).lower() in ("on", "off", "0", "1"):
            continue
        raw = Path(m.group(2))
        candidate = raw if raw.is_absolute() else cfg.parent / raw
        matches = []
        if candidate.is_file():
            files.append(candidate.resolve())
        elif candidate.parent.is_dir():
            matches = [
                p.resolve()
                for p in candidate.parent.glob(candidate.name + "*")
                if p.is_file()
            ]
            files.extend(matches)
        if not candidate.is_file() and not matches and missing is not None:
            missing.append(str(candidate))
    return files

def finite_number(row, key):
    try:
        return math.isfinite(float(row[key]))
    except (KeyError, TypeError, ValueError):
        return False

def inner_target_contract(row):
    rhs = float(row["rhs_rms"])
    target = float(row["target_absolute"])
    achieved = float(row["achieved_absolute"])
    if rhs <= 1e-300:
        return achieved <= target and achieved <= 1e-300
    return (math.isclose(target,
                         float(row["requested_relative_tolerance"]) * rhs,
                         rel_tol=5e-13, abs_tol=1e-300) and
            math.isclose(float(row["achieved_relative"]), achieved / rhs,
                         rel_tol=5e-13, abs_tol=1e-300))

def scientific_dirty_paths(status):
    critical_roots = ("lib/", "module/", "tools/", "CitcomS/Components/")
    critical_suffixes = (".c", ".h", ".py")
    paths = []
    for line in status.splitlines():
        path = line[3:].split(" -> ")[-1] if len(line) >= 4 else line
        if path.startswith(critical_roots) and path.endswith(critical_suffixes):
            paths.append(path)
    return paths

def manifest(args):
    repo = Path(args.repo).resolve()
    files = [Path(p).resolve() for p in args.input]
    missing_specs = []
    for cfg in [p for p in files if p.suffix == ".cfg" and p.is_file()]:
        files.extend(cfg_input_files(cfg, missing_specs))
    files = list(dict.fromkeys(files))
    binary = [Path(p).resolve() for p in args.binary]
    status = git(repo, "status", "--short")
    branch = git(repo, "symbolic-ref", "--short", "HEAD")
    head = git(repo, "rev-parse", "HEAD")
    mpicc = subprocess.getoutput("mpicc --version | head -1")
    config_status = repo / "config.status"
    build_command = os.environ.get("STAGE_ABC_BUILD_COMMAND", "").strip()
    if not build_command and config_status.is_file():
        try:
            build_command = subprocess.check_output(
                [str(config_status), "--config"], text=True,
                stderr=subprocess.STDOUT).strip()
        except (OSError, subprocess.CalledProcessError):
            build_command = ""
    compiler = os.environ.get("CC", "").strip() or mpicc
    scientific_dirty = scientific_dirty_paths(status)
    branch_ok = not args.expected_branch or branch == args.expected_branch
    provenance_complete = bool(head and compiler and build_command and branch_ok
                               and not scientific_dirty)
    result = {
        "schema": "strict-ala-stage-abc-v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source": {"repository": str(repo), "branch": branch,
                   "head": head, "expected_branch": args.expected_branch,
                   "branch_verified": branch_ok, "status": status,
                   "scientific_dirty": scientific_dirty,
                   "diff_binary": git(repo, "diff", "--binary")},
        "build": {"compiler": compiler,
                  "mpicc": mpicc,
                  "python": sys.version, "platform": platform.platform(),
                  "command": build_command,
                  "modules": os.environ.get("LOADEDMODULES", "")},
        "runtime": {k: os.environ.get(k, "") for k in
                    ("LSB_JOBID", "LSB_HOSTS", "OMP_NUM_THREADS", "MPI_LOCALRANKID")},
        "files": [{"path": str(p), "sha256": sha256(p)} for p in files + binary
                  if p.is_file()],
    }
    missing = [str(p) for p in files + binary if not p.is_file()] + missing_specs
    result["missing"] = missing
    result["provenance_complete"] = provenance_complete and not missing
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if missing or not result["provenance_complete"]:
        print("STAGE_A_FAIL missing=%s provenance_complete=%s scientific_dirty=%s" %
              (",".join(missing), result["provenance_complete"],
               ",".join(scientific_dirty)))
        return 1
    return 0

def generate_cfg(args):
    cfg_variant(args.base, args.c0, False)
    cfg_variant(args.base, args.c1, True)
    diff = normalized_cfg_diff(args.c0, args.c1)
    Path(args.diff).write_text("\n".join("%s: %s -> %s" % x for x in diff) + "\n")
    if [(x[0]) for x in diff] != ["ala_shallow_patch_preconditioner"]:
        print("STAGE_A_FAIL cfg_scientific_diff=" + repr(diff)); return 1
    return 0

def check_cfg(args):
    diff = normalized_cfg_diff(args.c0, args.c1)
    Path(args.diff).write_text("\n".join("%s: %s -> %s" % x for x in diff) + "\n")
    a, b = cfg_values(args.c0), cfg_values(args.c1)
    invalid_off = [k for k in REQUIRED_OFF if a.get(k) != "off" or b.get(k) != "off"]
    invalid_flags = (a.get("ala_stage_abc_production_logging") != "on" or
                     b.get("ala_stage_abc_production_logging") != "on" or
                     a.get("ala_stage_abc_adjoint_diagnostic") != "off" or
                     b.get("ala_stage_abc_adjoint_diagnostic") != "off" or
                     a.get("ala_shallow_patch_preconditioner") != "off" or
                     b.get("ala_shallow_patch_preconditioner") != "on")
    required_equal = ("ala_pressure_bpi_weight", "ala_augmented_lagrangian_gamma",
                      "ala_outer_solver", "ala_pcg_restart_interval",
                      "ala_inner_accuracy_max", "ala_inner_accuracy_factor",
                      "ala_unaugmented_momentum_tolerance", "tole_compressibility",
                      "nprocx", "nprocy", "nprocz", "nodex", "nodey", "nodez",
                      "levels", "steps")
    invalid_equal = [k for k in required_equal if not a.get(k) or a.get(k) != b.get(k)]
    required_values = {
        "ala_pressure_bpi_weight": "1.0",
        "ala_augmented_lagrangian_gamma": "10.0",
        "ala_outer_solver": "fgmres",
        "ala_pcg_restart_interval": "50",
        "ala_inner_accuracy_max": "1e-2",
        "ala_inner_accuracy_factor": "1e-2",
        "ala_unaugmented_momentum_tolerance": "1e-3",
        "tole_compressibility": "1e-02",
        "nprocx": "4", "nprocy": "4", "nprocz": "2",
        "nodex": "129", "nodey": "129", "nodez": "65",
        "levels": "5", "steps": "1",
    }
    invalid_values = [k for k, value in required_values.items()
                      if a.get(k) != value or b.get(k) != value]
    if ([x[0] for x in diff] != ["ala_shallow_patch_preconditioner"] or
            invalid_off or invalid_flags or invalid_equal or invalid_values):
        print("STAGE_A_FAIL cfg_scientific_diff=%r invalid_off=%r "
              "invalid_equal=%r invalid_values=%r" %
              (diff, invalid_off, invalid_equal, invalid_values))
        return 1
    return 0

def stage_b_decision(args):
    layouts = tuple(getattr(args, "expected_layout", None) or ("4x4x2",))
    with open(args.adjoint, newline="") as f: rows = list(csv.DictReader(f))
    expected_keys = {(layout, op, probe, velocity)
                     for layout in layouts for op in EXPECTED_OPERATORS
                     for probe in EXPECTED_PRESSURE_PROBES
                     for velocity in EXPECTED_VELOCITY_PROBES}
    keys = [(r.get("layout"), r.get("operator"), r.get("probe"),
             r.get("velocity_probe")) for r in rows]
    adjoint_complete = (set(keys) == expected_keys and len(keys) == len(set(keys))
                        and all(all(finite_number(r, k)
                                    for k in ADJOINT_NUMERIC_FIELDS) for r in rows))
    verdict = {}
    failure_keys = {}
    for op in ("D", "C", "B"):
        selected = [r for r in rows if r.get("operator") == op]
        failed = {(r["probe"], r["velocity_probe"], r["layout"])
                  for r in selected
                  if float(r["absolute_defect"]) > 100.0 * float(r["scale_floor"])
                  and float(r["delta_scale"]) > 1e-10}
        failure_keys[op] = sorted("%s/%s/%s" % x for x in failed)
        failed_modes = {(p, v) for p, v, layout in failed}
        common = {(p, v) for p, v in failed_modes
                  if all((p, v, layout) in failed for layout in layouts)}
        if not adjoint_complete: verdict[op] = "UNRESOLVED"
        elif failed and common != failed_modes: verdict[op] = "BOUNDARY_OR_OWNERSHIP_MISMATCH"
        elif failed: verdict[op] = "TRUE_DEFECT"
        elif any(float(r["denom_lr"]) <= float(r["scale_floor"]) for r in selected):
            verdict[op] = "NEAR_NULL_DENOMINATOR"
        else: verdict[op] = "PASS"
    with open(args.gauge, newline="") as f: gauge_rows = list(csv.DictReader(f))
    expected_gauge = {(layout, probe) for layout in layouts
                      for probe in EXPECTED_GAUGE_PROBES}
    gauge_keys = [(r.get("layout"), r.get("probe")) for r in gauge_rows]
    gauge_complete = (set(gauge_keys) == expected_gauge and
                      len(gauge_keys) == len(set(gauge_keys)) and
                      all(all(finite_number(r, k) for k in GAUGE_NUMERIC_FIELDS)
                          for r in gauge_rows))
    candidate_verdicts = {}
    for probe in ("P14_constant", "P15_density_gauge"):
        rr = [r for r in gauge_rows if r.get("probe") == probe]
        exact = gauge_complete and all(
            float(r["BTq_norm"]) <= 100.0 * float(r["scale_floor"]) and
            float(r["S_gamma_q_norm"]) <= 100.0 * float(r["scale_floor"])
            for r in rr)
        near = gauge_complete and all(
            float(r["ratio_B"]) <= 1.0e-2 and
            float(r["S_gamma_q_norm"]) <= 1.0e-3 * max(float(r["BTq_norm"]), 1e-300)
            for r in rr)
        candidate_verdicts[probe] = ("EXACT_NULL" if exact else
                                     "NEAR_NULL" if near else "NOT_NULL_CANDIDATE")
    if not gauge_complete: gauge = "UNRESOLVED"
    elif "EXACT_NULL" in candidate_verdicts.values(): gauge = "EXACT_NULL"
    elif "NEAR_NULL" in candidate_verdicts.values(): gauge = "NEAR_NULL"
    else: gauge = "UNRESOLVED"
    split_failed = [r for r in gauge_rows
                    if float(r["BT_split_absolute"]) > 100.0 * float(r["scale_floor"])
                    and float(r["BT_split_relative"]) > 1e-10] if gauge_complete else []
    if not gauge_complete:
        verdict["B"] = "UNRESOLVED"
    elif split_failed:
        verdict["B"] = "TRUE_DEFECT"
    blocked = any(verdict[x] in ("TRUE_DEFECT", "BOUNDARY_OR_OWNERSHIP_MISMATCH", "UNRESOLVED") for x in ("D", "C", "B"))
    out = {"D": verdict["D"], "C": verdict["C"], "B": verdict["B"],
           "GAUGE": gauge, "stage_C": "BLOCKED" if blocked else "ALLOWED",
           "gauge_candidates": candidate_verdicts,
           "adjoint_complete": adjoint_complete,
           "gauge_complete": gauge_complete,
           "expected_layouts": list(layouts),
           "failure_keys": failure_keys,
           "BT_split_failures": ["%s/%s" % (r["layout"], r["probe"])
                                 for r in split_failed],
           "schema": "strict-ala-stage-B-decision-v2"}
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print("STRICT_STAGE_B_BLOCKS_STAGE_C" if blocked else "STRICT_STAGE_B_ALLOWS_STAGE_C")
    return 1 if blocked else 0

def stage_c_decision(args):
    with open(args.iterations, newline="") as f:
        rows = list(csv.DictReader(f))
    with open(args.inner, newline="") as f:
        inner = list(csv.DictReader(f))
    out = {}
    valid = True
    terminal = {}
    iter_numeric = ("krylov_recursive", "krylov_explicit", "krylov_drift",
                    "continuity_numerator", "continuity_denominator",
                    "continuity_relative", "momentum_numerator",
                    "momentum_denominator", "momentum_relative",
                    "cumulative_inner_solves", "cumulative_inner_cycles",
                    "elapsed_seconds")
    for case in ("C0", "C1"):
        rr = [r for r in rows if r.get("case") == case]
        ir = [r for r in inner if r.get("case") == case]
        final_rows = [r for r in rr if r.get("final_iterate") == "1"]
        unique_iterations = len({r.get("iteration") for r in rr}) == len(rr)
        case_valid = (bool(rr) and bool(ir) and len(final_rows) == 1 and
                      unique_iterations and
                      all(all(finite_number(r, k) for k in iter_numeric) for r in rr))
        valid &= case_valid
        row = final_rows[0] if len(final_rows) == 1 else None
        viable = bool(case_valid and
                      float(row["continuity_relative"]) < 1e-2 and
                      float(row["momentum_relative"]) <= 1e-3)
        terminal[case] = row
        out[case] = {
            "viable": viable, "rows": len(rr), "inner_rows": len(ir),
            "last_iteration": int(row["iteration"]) if row else None,
            "continuity_relative": float(row["continuity_relative"]) if row else None,
            "momentum_relative": float(row["momentum_relative"]) if row else None,
            "inner_cycles": int(row["cumulative_inner_cycles"]) if row else None,
            "elapsed_seconds": float(row["elapsed_seconds"]) if row else None,
            "structure_valid": case_valid,
        }
    if not valid:
        decision = "INVALID_EXPERIMENT"
    elif out["C0"]["viable"] and not out["C1"]["viable"]:
        decision = "BPI_ONLY_WINS"
    elif out["C1"]["viable"] and not out["C0"]["viable"]:
        decision = "CONFIGURED_WINS"
    elif not out["C0"]["viable"] and not out["C1"]["viable"]:
        decision = "NO_PRODUCTION_APPROPRIATE_MAP"
    else:
        def dominates(a, b):
            return (a["inner_cycles"] <= 0.90 * b["inner_cycles"] and
                    a["elapsed_seconds"] <= 0.90 * b["elapsed_seconds"] and
                    a["continuity_relative"] <= 1.05 * b["continuity_relative"] and
                    a["momentum_relative"] <= 1.05 * b["momentum_relative"])
        if dominates(out["C0"], out["C1"]): decision = "BPI_ONLY_WINS"
        elif dominates(out["C1"], out["C0"]): decision = "CONFIGURED_WINS"
        else: decision = "TIE_NEEDS_REPEAT"
    out["decision"] = decision
    out["structure_valid"] = valid
    out["comparison_contract"] = {
        "minimum_savings_fraction": 0.10,
        "maximum_residual_regression_fraction": 0.05,
        "requires_both_wall_and_inner_cycle_savings": True,
    }
    out["schema"] = "strict-ala-stage-C-decision-v2"
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    return 0 if valid else 1

def final_audit(args):
    with open(args.iterations, newline="") as f:
        iterations = list(csv.DictReader(f))
    with open(args.inner, newline="") as f:
        inner = list(csv.DictReader(f))
    structural = True
    required_cases = {"C0", "C1"}
    try:
        iter_keys = [(r["case"], int(r["iteration"])) for r in iterations]
        structural &= (len(iter_keys) == len(set(iter_keys)) and bool(iter_keys)
                       and {r["case"] for r in iterations} == required_cases
                       and {r["case"] for r in inner} == required_cases)
        numeric_iter = ("krylov_recursive", "krylov_explicit", "krylov_drift",
                        "continuity_numerator", "continuity_denominator",
                        "continuity_relative", "momentum_numerator",
                        "momentum_denominator", "momentum_relative")
        structural &= all(all(finite_number(r, k) for k in numeric_iter)
                          for r in iterations)
        inner_numeric = ("rhs_rms", "requested_relative_tolerance",
                         "target_absolute", "achieved_absolute",
                         "achieved_relative", "seconds")
        structural &= bool(inner) and all(
            r["status"] == "CONVERGED" and float(r["achieved_absolute"]) <= float(r["target_absolute"])
            and int(r["cycles"]) < int(r["max_cycles"])
            and int(r["cycles"]) >= 0
            and all(finite_number(r, k) and float(r[k]) >= 0.0 for k in inner_numeric)
            for r in inner)
        for case in required_cases:
            rr = [r for r in iterations if r["case"] == case]
            ir = [r for r in inner if r["case"] == case]
            call_ids = sorted(int(r["call_id"]) for r in ir)
            structural &= (sum(r.get("final_iterate") == "1" for r in rr) == 1
                           and int(rr[-1]["cumulative_inner_solves"]) == len(ir)
                           and call_ids == list(range(1, len(ir)+1))
                           and all(inner_target_contract(r) for r in ir))
    except (KeyError, TypeError, ValueError): structural = False
    with open(args.manifest) as f:
        manifest_data = json.load(f)
    input_identity = all(Path(x["path"]).is_file() and sha256(x["path"]) == x["sha256"]
                         for x in manifest_data.get("files", []))
    provenance_complete = bool(manifest_data.get("provenance_complete"))
    structural &= input_identity and not manifest_data.get("missing") and provenance_complete
    with open(args.stage_b_decision) as f:
        stage_b = json.load(f)
    stage_b_pass = (stage_b.get("stage_C") == "ALLOWED" and
                    all(stage_b.get(k) == "PASS" for k in ("D", "C", "B")) and
                    stage_b.get("adjoint_complete") is True and
                    stage_b.get("gauge_complete") is True)
    cfg_diff_verified = Path(args.cfg_diff).read_text().strip() == \
        "ala_shallow_patch_preconditioner: off -> on"
    binary_unchanged = Path(args.binary_pre).read_bytes() == Path(args.binary_post).read_bytes()
    hidden_fallback = False
    completed_log_cases = set()
    for log in args.case_log:
        text = Path(log).read_text(errors="replace")
        hidden_fallback |= bool(re.search(r"(?:fallback_blocks|velocity_block_fallbacks)=[1-9][0-9]*", text))
        completed_log_cases.update(re.findall(
            r"STRICT_ALA_STAGE_C_CASE_COMPLETE case=(C0|C1) ", text))
    case_logs_complete = (len(args.case_log) == 2 and
                          completed_log_cases == required_cases)
    structural &= (stage_b_pass and cfg_diff_verified and binary_unchanged and
                   not hidden_fallback and case_logs_complete)
    costs = []
    for case in ("C0", "C1"):
        rr = [r for r in inner if r.get("case") == case]
        ir = [r for r in iterations if r.get("case") == case]
        costs.append({"case": case, "K_gamma_inverse_calls": len(rr),
                      "K_gamma_applications": int(ir[-1]["cumulative_K_gamma_applications"]) if ir else 0,
                      "schur_actions": int(ir[-1]["cumulative_schur_actions"]) if ir else 0,
                      "preconditioner_applications": int(ir[-1]["cumulative_preconditioner_applications"]) if ir else 0,
                      "bpi_construction_seconds": float(ir[-1]["bpi_construction_seconds"]) if ir else 0.0,
                      "schwarz_construction_seconds": float(ir[-1]["schwarz_construction_seconds"]) if ir else 0.0,
                      "inner_cycles": sum(int(r["cycles"]) for r in rr),
                      "inner_seconds": sum(float(r["seconds"]) for r in rr),
                      "fgmres_seconds": float(ir[-1]["elapsed_seconds"]) if ir else 0.0})
    with open(args.cost, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=costs[0].keys()); w.writeheader(); w.writerows(costs)
    with open(args.decision) as f:
        decision = json.load(f)
    structural &= (decision.get("structure_valid") is True and
                   decision.get("decision") != "INVALID_EXPERIMENT")
    final = {"stage_A": "PASS" if provenance_complete and input_identity else "FAIL",
             "stage_B": "PASS" if stage_b_pass else "FAIL",
             "stage_C0": "PASS_PHYSICAL" if decision["C0"]["viable"] else "NONCONVERGED",
             "stage_C1": "PASS_PHYSICAL" if decision["C1"]["viable"] else "NONCONVERGED",
             "decision": decision["decision"] if structural else "INVALID_EXPERIMENT",
             "provenance_complete": provenance_complete,
             "binary_unchanged": binary_unchanged,
             "input_identity_verified": input_identity,
             "cfg_single_variable_verified": cfg_diff_verified,
             "case_logs_complete": case_logs_complete,
             "hidden_fallback_detected": hidden_fallback,
             "structural_checks_pass": structural,
             "schema": "strict-ala-stage-ABC-final-v2"}
    Path(args.output).write_text(json.dumps(final, indent=2, sort_keys=True) + "\n")
    return 0 if structural else 1

def main():
    p = argparse.ArgumentParser(); s = p.add_subparsers(dest="cmd", required=True)
    g = s.add_parser("generate-cfg"); g.add_argument("--base", required=True); g.add_argument("--c0", required=True); g.add_argument("--c1", required=True); g.add_argument("--diff", required=True); g.set_defaults(fn=generate_cfg)
    k = s.add_parser("check-cfg"); k.add_argument("--c0", required=True); k.add_argument("--c1", required=True); k.add_argument("--diff", required=True); k.set_defaults(fn=check_cfg)
    m = s.add_parser("manifest"); m.add_argument("--repo", required=True); m.add_argument("--expected-branch", default=""); m.add_argument("--output", required=True); m.add_argument("--input", nargs="*", default=[]); m.add_argument("--binary", nargs="*", default=[]); m.set_defaults(fn=manifest)
    b = s.add_parser("decide-stage-b"); b.add_argument("--adjoint", required=True); b.add_argument("--gauge", required=True); b.add_argument("--expected-layout", action="append"); b.add_argument("--output", required=True); b.set_defaults(fn=stage_b_decision)
    c = s.add_parser("decide-stage-c"); c.add_argument("--iterations", required=True); c.add_argument("--inner", required=True); c.add_argument("--output", required=True); c.set_defaults(fn=stage_c_decision)
    f = s.add_parser("final-audit"); f.add_argument("--iterations", required=True); f.add_argument("--inner", required=True); f.add_argument("--cost", required=True); f.add_argument("--decision", required=True); f.add_argument("--manifest", required=True); f.add_argument("--stage-b-decision", required=True); f.add_argument("--cfg-diff", required=True); f.add_argument("--binary-pre", required=True); f.add_argument("--binary-post", required=True); f.add_argument("--case-log", action="append", required=True); f.add_argument("--output", required=True); f.set_defaults(fn=final_audit)
    a = p.parse_args(); return a.fn(a)
if __name__ == "__main__": sys.exit(main())
