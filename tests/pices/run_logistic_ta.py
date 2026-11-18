"""Small MPI integration gate for depth-evaluated logistic particle TA.

Uses PICESManufactured, the real production kernels with prescribed zero flow.
Runs only under the explicitly supplied output directory (normally /tmp).
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess

import numpy as np


CENTER = 80.0 / 127.42
WIDTH = 4.0 / 127.42
RANKS = 12


def fields(line):
    return dict(word.split("=", 1) for word in line.split() if "=" in word)


def replace_settings(base, changes):
    for key, value in changes.items():
        pattern = r"^" + re.escape(key) + r"\s*=.*$"
        line = key + "=" + str(value)
        base = (re.sub(pattern, lambda match: line, base, flags=re.M)
                if re.search(pattern, base, re.M) else base + "\n" + line + "\n")
    return base


def local_thickness(flag, maximum):
    distance = np.abs(flag)
    factor = np.where(distance < .06, .003,
                      np.where(distance <= .1, .003 + .997 * (distance - .06) / .04, 1.))
    signed_side = (flag * 1000 / 6371 >= -.1) & (flag <= 0) & (distance * 1000 / 6371 <= .1)
    return maximum * np.where(signed_side, .003, factor)


def quadrature_rate(fraction, center, width):
    """Independent quadrature of the defining integral, not the C antiderivative."""
    knots, weights = np.polynomial.legendre.leggauss(128)
    bottom = 1 / (1 + math.exp((1 - center) / width))
    def integral(left):
        points = left[..., None] + (1 - left[..., None]) * (knots + 1) / 2
        values = 1 / (1 + np.exp((points - center) / width)) - bottom
        return (1 - left) / 2 * (values @ weights)
    return integral(fraction) / integral(np.asarray(0.))


def accepted_records(directory):
    lines = (directory / "DATA/0/log").read_text().splitlines()
    diagnostics = [fields(line) for line in lines if line.startswith("PICES_DIAGNOSTICS ")]
    ta = [fields(line) for line in lines if line.startswith("PICES_TA ")]
    return lines, diagnostics, ta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path)
    parser.add_argument("cfg", type=Path, help="standalone P4 fixture")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    executable = args.build.resolve() / "PICESManufactured"
    base = args.cfg.read_text()
    common = dict(nodex=5, nodey=5, nodez=5, mgunitx=4, mgunity=4, mgunitz=4,
                  levels=1, nprocx=1, nprocy=1, nprocz=1, maxstep=1, maxtotstep=2,
                  storage_spacing=1, checkpointFrequency=1, tracers_per_element=32,
                  pices_checkpoint="on", pices_projection="lumped", Q0=0,
                  dissipation_number=0, phase_delta_rho="0,0,0", phase_delta_s="0,0,0",
                  fixed_timestep="3.77e-8", kT_exponent=0, kd_upper_prefactor=4,
                  kd_lower_prefactor=4, lith_age=1, lith_age_asml=1,
                  lith_age_asml_profile=1, lith_age_asml_logistic_center=CENTER,
                  lith_age_asml_logistic_width=WIDTH, lith_age_depth=.02,
                  lith_age_asml_tau_Ma=.25, max_plate_age_Ma=70, start_age=2.13,
                  pices_ta_before_heat="off", restart="off", tracer_seed_distribution="legacy")

    def run(name, changes=None, expected_failure=False, invalid_input=False):
        directory = output / name
        directory.mkdir()
        (directory / "DATA").mkdir()
        config = replace_settings(base, dict(common, **(changes or {})))
        (directory / "case.cfg").write_text(config)
        reference = re.search(r"^refstate_file\s*=\s*(.+)$", config, re.M)[1]
        (directory / reference).write_text("1 1 .6 1 1\n" * 5)
        forcing = directory / "pices_p4_forcing"
        forcing.mkdir()
        for age in range(4):
            (forcing / ("trench.%d.xyz" % age)).write_text("")
            for cap in range(12):
                (forcing / ("age.%d.%d" % (age, cap))).write_text("40\n" * 25)
        with (directory / "stdout").open("w") as stdout, (directory / "stderr").open("w") as stderr:
            result = subprocess.run(
                [os.environ.get("MPIEXEC", "/usr/local/bin/mpiexec"), "--oversubscribe", "-n", str(RANKS),
                 str(executable), "case.cfg"], cwd=directory,
                env=dict(os.environ, PICES_TEST_MODE="constant"), stdout=stdout, stderr=stderr, timeout=180)
        error = (directory / "stderr").read_text()
        if invalid_input:
            assert result.returncode != 0 and "lith_age_asml_profile must be" in error, (name, result.returncode, error[-2000:])
            assert not any("PICES_STEP " in path.read_text() for path in (directory / "DATA").glob("*/log"))
        elif expected_failure:
            assert result.returncode == 72 and "metadata/physics" in error, (name, result.returncode, error[-2000:])
            assert not any("PICES_STEP " in path.read_text() for path in (directory / "DATA").glob("*/log"))
        else:
            assert result.returncode in (0, 8), (name, result.returncode, error[-2000:])
            assert "PICES_ERROR" not in error, (name, error[-2000:])
            final_step = int((changes or {}).get("maxstep", common["maxstep"]))
            assert all((directory / ("particles.%d.%d.txt" % (rank, final_step))).exists()
                       for rank in range(RANKS)), name
        print(name + ": PASS", flush=True)
        return directory

    off = run("after_heat_off", dict(lith_age=0, lith_age_asml=0))
    on = run("after_heat_on")
    maximum_error = 0.
    supported_count = changed_count = below_count = nodal_active = 0
    ta_energy = {}
    grid_energy = {}
    center, width = float(np.float32(CENTER)), float(np.float32(WIDTH))
    for rank in range(RANKS):
        old = np.loadtxt(off / ("particles.%d.1.txt" % rank))
        new = np.loadtxt(on / ("particles.%d.1.txt" % rank))
        assert np.array_equal(old[:, 1:], new[:, 1:]), "TA changed positions/weights"
        nodes_old = np.loadtxt(off / ("nodes.%d.1.txt" % rank))
        nodes_new = np.loadtxt(on / ("nodes.%d.1.txt" % rank))
        path = on / ("ta.%d.1.txt" % rank)
        header = {key: float(value) for key, value in fields(path.read_text().splitlines()[0]).items()}
        data = np.loadtxt(path)
        ids = new[:, 1:17:2].astype(int) - 1
        weights = new[:, 2:17:2]
        depth = 1 - np.linalg.norm(new[:, 17:20], axis=1)
        thickness = local_thickness(data[:, 5], header["depth"])
        nodal_active += int(np.count_nonzero((data[:, 2] < thickness) & (data[:, 2] > 0) & (nodes_new[:, 5] == 0)))
        h = thickness[ids]
        active = (depth[:, None] > 0) & (depth[:, None] < h)
        fraction = np.zeros_like(h)
        physical_dt = float(np.float32(np.float32(header["dt"]) * np.float32(header["scalet"])))
        x = np.broadcast_to(depth[:, None], h.shape)[active] / h[active]
        fraction[active] = -np.expm1(-physical_dt * quadrature_rate(x, center, width) / header["tau"])
        capped_age = np.minimum(data[:, 4], float(np.float32(header["cap"] / header["scalet"])))
        eta = np.divide(.5 * depth[:, None], np.sqrt(np.maximum(capped_age[ids], 1.e-300)))
        eta = np.where(capped_age[ids] > 0, eta, 10.)
        erfc = np.fromiter((math.erfc(value) for value in eta.flat), float, eta.size).reshape(eta.shape)
        targets = data[ids, 3] - header["surface"] * erfc
        reaction = weights * fraction
        alpha = reaction.sum(axis=1)
        expected = (1 - alpha) * old[:, 0] + (reaction * targets).sum(axis=1)
        maximum_error = max(maximum_error, float(np.max(np.abs(new[:, 0] - expected))))
        supported = alpha > 0
        supported_count += int(supported.sum())
        changed_count += int(np.count_nonzero(new[supported, 0] != old[supported, 0]))
        below = depth >= header["depth"]
        below_count += int(below.sum())
        assert np.array_equal(new[below, 0], old[below, 0]), "TA leaked below H"
        assert np.all(new[:, 0] >= np.minimum(old[:, 0], targets.min(axis=1)) - 2.e-14)
        assert np.all(new[:, 0] <= np.maximum(old[:, 0], targets.max(axis=1)) + 2.e-14)
        for xyz, mass, delta, total_delta in zip(nodes_new[:, 1:4], nodes_new[:, 6], data[:, 1],
                                                nodes_new[:, 4] - nodes_old[:, 4]):
            key = tuple(np.round(xyz, 10))
            ta_energy[key] = mass * delta
            grid_energy[key] = mass * total_delta
        for side in ("surf", "botm"):
            relative = "DATA/%d/q.%s.%d.1" % (rank, side, rank)
            assert (on / relative).read_bytes() == (off / relative).read_bytes(), "after-heat TA changed CBF"
    assert nodal_active == 0 and supported_count > 0 and changed_count > 0 and below_count > 0
    assert maximum_error < 3.e-14, maximum_error
    _, diagnostics, ta = accepted_records(on)
    assert len(ta) == 1 and ta[0]["calls"] == "1" and ta[0]["stage"] == "after_heat"
    assert ta[0]["operator"] == "particle_relaxation_logistic_v1"
    closure = max(abs(float(row["closure_residual_nd"])) for row in diagnostics)
    assert closure < 1.e-12, closure
    storage_error = abs(sum(ta_energy.values()) - float(ta[0]["storage"]))
    grid_error = abs(sum(grid_energy.values()) - float(ta[0]["total_grid_storage"]))
    assert max(storage_error, grid_error) < 1.e-12, (storage_error, grid_error)

    preheat = dict(pices_ta_before_heat="on", maxstep=2, maxtotstep=3)
    continuous = run("continuous", preheat)
    split = run("split", dict(pices_ta_before_heat="on"))
    resume = dict(preheat, restart="on", solution_cycles_init=1,
                  datadir_old=str(split / "DATA/%RANK"), datafile_old="PICES_P4")
    restarted = run("restart", resume)
    for directory, expected_steps in ((continuous, [1, 2]), (split, [1]), (restarted, [2])):
        _, records, events = accepted_records(directory)
        assert [int(event["step"]) for event in events] == expected_steps
        assert all(event["calls"] == "1" and event["stage"] == "before_heat" for event in events)
        assert all(abs(float(row["closure_residual_nd"])) < 1.e-12 for row in records)
    for rank in range(RANKS):
        for stem in ("particles", "nodes"):
            relative = "%s.%d.2.txt" % (stem, rank)
            assert (continuous / relative).read_bytes() == (restarted / relative).read_bytes(), relative
        relative = "DATA/%d/PICES_P4.chkpt.%d.2.pices.state" % (rank, rank)
        assert (continuous / relative).read_bytes() == (restarted / relative).read_bytes(), relative
    rejected = {}
    for name, key, value in (("profile", "lith_age_asml_profile", 0),
                             ("center", "lith_age_asml_logistic_center", .6),
                             ("width", "lith_age_asml_logistic_width", .04),
                             ("tau", "lith_age_asml_tau_Ma", .3)):
        run("reject_" + name, dict(resume, **{key: value}), expected_failure=True)
        rejected[name] = "PASS"
    invalid = {}
    for name, key, value in (("profile", "lith_age_asml_profile", 2),
                             ("center", "lith_age_asml_logistic_center", 0),
                             ("width", "lith_age_asml_logistic_width", 0)):
        run("invalid_" + name, {key: value}, invalid_input=True)
        invalid[name] = "PASS"
    summary = dict(status="PASS", ranks=RANKS, particle_reaction_max_error_nd=maximum_error,
                   supported_particles=supported_count, changed_particles=changed_count,
                   particles_below_H_unchanged=below_count, eligible_grid_nodes=nodal_active,
                   closure_max_nd=closure, TA_storage_error_nd=storage_error,
                   grid_storage_error_nd=grid_error, CBF_after_heat_invariant=True,
                   preheat_two_step_restart_byte_identical=True, restart_rejections=rejected,
                   invalid_input_rejections=invalid)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
