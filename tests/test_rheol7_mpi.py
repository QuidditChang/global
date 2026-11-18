"""Run each reference rheology through one real 12-rank PICES/pPv step.

Build with tests/cbf/build_validation.py. Supply the small standalone
cmbhf_EBA_PICES_P3.cfg fixture, not a production Pyre configuration.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("cfg", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    base = {}
    for line in args.cfg.read_text().splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            base[key.strip()] = value.strip()
    common = dict(maxstep=1, maxtotstep=2, checkpointFrequency=1,
                  tracers_per_element=32, pices_checkpoint="on", pices_p4="on",
                  Q0=0, dissipation_number=0, phase_delta_rho="0,0,0",
                  phase_delta_s="0,0,0", fixed_timestep="1e-7",
                  cold_scale=1, VMIN="on", VMAX="on", visc_min=.1, visc_max=100,
                  refvisc="1e21", rheol7_layer_factors="2,.25,.5,1.5",
                  pPv_enabled="on", pPv_depth=.25, pPv_transT=.6,
                  pPv_clapeyron=.1, pPv_width=.03, pPv_delta_rho=66.4,
                  pPv_delta_s=-.014, pPv_latent_enabled="off",
                  pPv_k_ratio=1, pPv_eta_ratio=1)
    report = {}
    for mode in (7, 71, 72):
        run = output / str(mode)
        run.mkdir()
        (run / "DATA").mkdir()
        values = dict(base, **common)
        values["rheol"] = mode
        (run / "case.cfg").write_text("".join("%s=%s\n" % item for item in values.items()))
        (run / values["refstate_file"]).write_text("1 1 .6 1 1\n" * int(values["nodez"]))
        with (run / "stdout").open("w") as stdout, (run / "stderr").open("w") as stderr:
            result = subprocess.run(
                [os.environ.get("MPIEXEC", "/usr/local/bin/mpiexec"), "--oversubscribe",
                 "-n", "12", str(args.build.resolve() / "CitcomSFull"), "case.cfg"],
                cwd=run, stdout=stdout, stderr=stderr, timeout=180)
        errors = (run / "stderr").read_text()
        assert result.returncode in (0, 8), (mode, result.returncode, errors[-3000:])
        assert "PICES_ERROR" not in errors and "rheol_family_v1" in errors, errors[-3000:]
        assert "fatal:" not in errors, errors[-3000:]
        assert "Average cpu time taken for velocity step" in errors, errors[-3000:]
        assert "rheol=%d: TDEPV=1 cold_scale=1" % mode in errors
        for rank in range(12):
            log = (run / ("DATA/%d/log" % rank)).read_text()
            assert "PICES_STEP step=1 " in log, (mode, rank, log[-2000:])
        metadata = list((run / "DATA").glob("*/*.pices.json"))
        assert len(metadata) >= 12, (mode, "missing checkpoint metadata", list((run / "DATA/0").iterdir()))
        report[str(mode)] = {"status": "PASS", "ranks": 12, "accepted_steps": 1,
                             "pPv_enabled": True, "checkpoint_sidecars": len(metadata)}
        print("PASS: rheol=%d, 12 ranks, coupled PICES/pPv step and checkpoint" % mode, flush=True)
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
