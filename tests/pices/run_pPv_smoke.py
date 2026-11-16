"""12-rank pPv integration: prescribed heat comparisons and a coupled step.

Build with tests/cbf/build_validation.py and tests/pices/build_manufactured.py.
The cfg argument is the small standalone cmbhf_EBA_PICES_P3.cfg fixture, not
the production Pyre configuration. All test inputs/outputs stay in --output.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("cfg", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    base = {}
    for line in args.cfg.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            base[key.strip()] = value.strip()
    common = dict(maxstep=1, maxtotstep=2, tracers_per_element=32,
                  pices_checkpoint="off", Q0=0, dissipation_number=0,
                  phase_delta_rho="0,0,0", phase_delta_s="0,0,0",
                  fixed_timestep="1e-7", profile_optional="temperature,density_total,k,qphase",
                  pPv_depth=.25, pPv_transT=.6, pPv_clapeyron=.1,
                  pPv_width=.03, pPv_delta_rho=66.4, pPv_delta_s=-.014,
                  pPv_latent_enabled="off", pPv_k_ratio=1, pPv_eta_ratio=1)
    results = {}

    def run(name, changes, coupled=False):
        directory = root / name
        directory.mkdir()
        (directory / "DATA").mkdir()
        values = {**base, **common, **changes}
        (directory / "case.cfg").write_text("".join(f"{k}={v}\n" for k, v in values.items()))
        (directory / values["refstate_file"]).write_text("1 1 .6 1 1\n" * int(values["nodez"]))
        executable = "CitcomSFull" if coupled else "PICESManufactured"
        with (directory / "stdout").open("w") as out, (directory / "stderr").open("w") as err:
            result = subprocess.run(
                [os.environ.get("MPIEXEC", "mpiexec"), "--oversubscribe", "-n", "12",
                 str(args.build.resolve() / executable), "case.cfg"],
                cwd=directory, stdout=out, stderr=err, timeout=120,
                env=dict(os.environ, PICES_TEST_MODE="diffusion", PICES_TEST_DIAGNOSTICS="1"))
        # Legacy normal finalization exits with status 8.
        assert result.returncode in (0, 8), (name, result.returncode,
                                              (directory / "stderr").read_text()[-3000:])
        for rank in range(12):
            log = (directory / f"DATA/{rank}/log").read_text()
            assert "PICES_STEP step=1 " in log, (name, rank)
        if not coupled:
            data = np.concatenate([np.loadtxt(directory / f"nodes.{rank}.1.txt")[:, 4]
                                   for rank in range(12)])
            assert np.all(np.isfinite(data))
        else:
            data = None
        results[name] = {"status": "PASS", "coupled_stokes": coupled}
        print(name, "PASS", flush=True)
        return directory, data

    off_dir, off = run("off", {"pPv_enabled": "off"})
    density_dir, density = run("density", {"pPv_enabled": "on"})
    latent_dir, latent = run("latent", {"pPv_enabled": "on", "pPv_latent_enabled": "on"})
    transport_dir, transport = run("transport", {"pPv_enabled": "on", "pPv_k_ratio": 1.5,
                                                 "pPv_eta_ratio": .1})
    assert np.array_equal(off, density), "density-only changed prescribed-flow heat evolution"
    assert np.max(abs(latent - density)) > 1e-12, "latent switch did not reach heat evolution"
    assert np.max(abs(transport - density)) > 1e-12, "conductivity ratio did not reach heat evolution"
    for name, directory, enabled in (("off", off_dir, False), ("density", density_dir, True)):
        profiles = sorted(directory.rglob("*.npz"))
        assert profiles, (name, "profile output absent")
        keys = set()
        for path in profiles:
            with np.load(path) as profile:
                keys.update(profile.files)
        assert any("pPv_fraction" in key for key in keys) == enabled, (name, keys)
    results["prescribed_comparison"] = {
        "density_only_exactly_equal": True,
        "latent_max_temperature_change_nd": float(np.max(abs(latent-density))),
        "conductivity_max_temperature_change_nd": float(np.max(abs(transport-density))),
    }
    run("coupled", {"pPv_enabled": "on", "pPv_latent_enabled": "on",
                    "pPv_k_ratio": 1.5, "pPv_eta_ratio": .1}, coupled=True)
    (root / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
