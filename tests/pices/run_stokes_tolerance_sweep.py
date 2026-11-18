"""Exercise configurable Stokes tolerances on the small 12-rank MG fixture.

This checks parameter plumbing, actual targets, initialization identity and
equivalence of implicit and explicit D; it is not a production performance
benchmark.
"""
import argparse
import gzip
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess


def fields(line):
    return dict(re.findall(r"(\w+)=([^\s]+)", line))


def edit(text, changes):
    for key, value in changes.items():
        pattern = r"^" + re.escape(key) + r"\s*=.*$"
        line = "%s=%s" % (key, value)
        if re.search(pattern, text, re.M):
            text = re.sub(pattern, lambda _: line, text, flags=re.M)
        else:
            text += "\n" + line + "\n"
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("runs", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = args.runs.resolve() / "pices_p8b"
    base = edit((source / "multigrid.cfg").read_text(), {
        "maxstep": 1, "maxtotstep": 2, "pices_checkpoint": "off",
        "CBF_frequency": 0, "output_q_surf_CBF": "off", "output_q_botm_CBF": "off",
        "pPv_enabled": "on", "pPv_depth": .25, "pPv_transT": .6,
        "pPv_clapeyron": .1, "pPv_width": .03, "pPv_delta_rho": 66.4,
        "pPv_delta_s": -.014, "pPv_latent_enabled": "off",
        "pPv_k_ratio": 1, "pPv_eta_ratio": 1,
    })
    summary = {}

    def run(name, changes, invalid=None):
        directory = root / name
        directory.mkdir()
        (directory / "DATA").mkdir()
        (directory / "case.cfg").write_text(edit(base, changes))
        shutil.copy2(source / "refstate.txt", directory)
        shutil.copytree(source / "forcing", directory / "forcing",
                        ignore=shutil.ignore_patterns("._*"))
        env = {k: v for k, v in os.environ.items() if not k.startswith("PICES_")}
        with (directory / "stdout").open("w") as out, (directory / "stderr").open("w") as err:
            proc = subprocess.run([
                os.environ.get("MPIEXEC", "/usr/local/bin/mpiexec"), "--oversubscribe",
                "-n", "12", str(args.build.resolve() / "CitcomSFull"), "case.cfg",
            ], cwd=directory, env=env, stdout=out, stderr=err, timeout=180)
        log_path = directory / "DATA/0/log"
        log = log_path.read_text() if log_path.exists() else ""
        stderr = (directory / "stderr").read_text()
        if invalid:
            assert proc.returncode != 0 and invalid in stderr, (name, proc.returncode, stderr[-2000:])
            assert "PICES_STOKES_INIT" not in log
            summary[name] = {"status": "PASS", "rejected": invalid}
            print(name, "rejected", flush=True)
            return directory, []
        policy_failure = (name in ("B", "C", "D", "E") and proc.returncode == 72
                          and "P3 " in stderr and "did not converge" in stderr)
        assert proc.returncode in (0, 8) or policy_failure, (name, proc.returncode, stderr[-3000:])
        stokes = [fields(l) for l in log.splitlines() if l.startswith("PICES_STOKES ")]
        expected_policy = changes.get("pices_stokes_inner_policy", "legacy_initial_force")
        assert all(r["policy"] == expected_policy for r in stokes), (name, stokes)
        if not policy_failure:
            assert len(stokes) >= 2 and all(r["status"] == "PASS" for r in stokes), (name, stokes)
        inner = [fields(l) for l in log.splitlines() if l.startswith("PICES_STOKES_INNER ")]
        identity = [fields(l) for l in log.splitlines() if l.startswith("PICES_STOKES_INIT ")]
        diagnostic = changes.get("pices_stokes_diagnostics") == "on"
        assert bool(inner) == diagnostic and bool(identity) == diagnostic
        for row in inner:
            assert row["policy"] == expected_policy, (name, row)
            if row["status"] == "PASS":
                assert row["actual_status"] == "PASS", (name, row)
                assert float(row["actual"]) < float(row["target"])
            else:
                assert policy_failure, (name, row)
            if row["policy"] == "current_rhs":
                expected = max(1e-14, float(changes["pices_stokes_inner_rtol"]) * float(row["rhs_rms"]))
            else:
                outer = next((r for r in stokes if r["step"] == row["step"]), None)
                if outer is None:
                    # A failing initial correction has no outer summary yet;
                    # its existing force line has nine decimal digits.
                    force = re.findall(r"initial residue of momentum equation F ([^ ]+)", log)
                    expected = 1e-6 * float(force[int(row["step"])])
                    assert math.isclose(float(row["target"]), expected, rel_tol=2e-9)
                    continue
                expected = float(outer["outer_accuracy"]) * float(outer["initial_force_rms"])
            assert math.isclose(float(row["target"]), expected, rel_tol=2e-14, abs_tol=1e-25), (name, row, expected)
        summary[name] = {"status": "PASS", "policy_result": "EXPECTED_POLICY_FAILURE" if policy_failure else "CONVERGED",
                         "inner_solves": len(inner),
                         "outer_iterations": [int(r["iterations"]) for r in stokes],
                         "initial_identity": identity[0] if identity else None}
        print(name, "PASS", flush=True)
        return directory, identity

    default, _ = run("default", {})
    explicit, _ = run("explicit_default", {"pices_stokes_inner_policy": "legacy_initial_force",
                                            "pices_stokes_inner_rtol": "1e-6"})
    compared = 0
    for path in default.glob("DATA/*/*/*.gz"):
        other = explicit / path.relative_to(default)
        assert other.exists()
        assert gzip.open(path, "rb").read() == gzip.open(other, "rb").read(), path
        compared += 1
    assert compared > 0
    identities = []
    for name, policy, rtol in (("A", "current_rhs", "1e-6"),
                               ("B", "current_rhs", "1e-5"),
                               ("C", "current_rhs", "1e-4"),
                               ("D", "legacy_initial_force", "1e-6"),
                               ("E", "current_rhs", "1e-7")):
        directory, identity = run(name, {"pices_stokes_inner_policy": policy,
                                        "pices_stokes_inner_rtol": rtol,
                                        "pices_stokes_diagnostics": "on"})
        identities.append(identity[0])
        if name == "D":
            for path in default.glob("DATA/*/*/*.gz"):
                assert gzip.open(path, "rb").read() == gzip.open(directory / path.relative_to(default), "rb").read(), path
    assert all(row == identities[0] for row in identities), identities
    for name, changes, message in (
        ("bad_policy", {"pices_stokes_inner_policy": "unknown"}, "pices_stokes_inner_policy must be"),
        ("long_policy", {"pices_stokes_inner_policy": "x" * 300}, "pices_stokes_inner_policy is too long"),
        ("zero_rtol", {"pices_stokes_inner_rtol": "0"}, "pices_stokes_inner_rtol must be"),
        ("nan_rtol", {"pices_stokes_inner_rtol": "nan"}, "pices_stokes_inner_rtol must be"),
    ):
        run(name, changes, invalid=message)
    summary["default_regression"] = {"status": "PASS", "decoded_files_equal": compared,
                                      "default_policy": "legacy_initial_force",
                                      "diagnostics_do_not_change_outputs": True,
                                      "initial_identity_equal_all_policies": True}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
