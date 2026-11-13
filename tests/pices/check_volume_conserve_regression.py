"""Compare the default-off solver with an independent archived HEAD baseline.

Only the temporary baseline build helper's commit lookup is substituted. Solver
sources are exactly git-archive HEAD, with no source edits and no production
build-file changes. Each binary also passes the existing real off/on audit test.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[2]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def run(command, logfile):
    with logfile.open('w') as stream:
        completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError('Command failed; inspect ' + str(logfile))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('current_build', type=Path)
    parser.add_argument('runs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--radial-partitions', type=int, choices=(1, 2), default=1)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    commit = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    archive = output / 'baseline.tar'
    subprocess.run(['git', '-C', str(ROOT), 'archive', '--format=tar', '-o', str(archive),
                    commit, 'lib', 'bin', 'tests/cbf'], check=True)
    source = output / 'baseline_source'
    source.mkdir()
    with tarfile.open(archive) as bundle:
        for member in bundle.getmembers():
            target = source / member.name
            if not target.resolve().is_relative_to(source) or member.issym() or member.islnk():
                raise ValueError('Unexpected path/link in source archive')
        bundle.extractall(source)
    helper = source / 'tests/cbf/build_validation.py'
    text = helper.read_text()
    lookup = "subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip()"
    if text.count(lookup) != 1:
        raise ValueError('Baseline helper git lookup changed; do not apply an unreviewed rewrite')
    helper.write_text(text.replace(lookup, repr(commit)))
    baseline_build = output / 'baseline_build'
    frozen_build = output / 'current_build'
    frozen_build.mkdir()
    current_binary = args.current_build.resolve() / 'CitcomSFull'
    original_hash = digest(current_binary)
    shutil.copy2(current_binary, frozen_build / 'CitcomSFull')
    if digest(frozen_build / 'CitcomSFull') != original_hash:
        raise ValueError('Current solver changed while freezing binary')
    run([sys.executable, str(helper), '--build-dir', str(baseline_build)], output / 'baseline_build.log')
    test = ROOT / 'tests/pices/check_volume_audit_integration.py'
    tests = {}
    for version, build in [('baseline', baseline_build), ('current', frozen_build)]:
        run([sys.executable, str(test), str(build), str(args.runs.resolve()), str(output / version),
             '--radial-partitions', str(args.radial_partitions)], output / (version + '_audit.log'))
        tests[version] = json.loads((output / version / 'results.json').read_text())
        if tests[version]['status'] != 'PASS':
            raise ValueError('Audit integration did not pass: ' + version)
    spec = importlib.util.spec_from_file_location('volume_audit_integration', test)
    integration = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(integration)
    ranks = 12 * args.radial_partitions
    compared = 0
    for audit_mode in ('off', 'on'):
        for rank in range(ranks):
            for step in (0, 2):
                for suffix in ('', '.pices.state'):
                    relative = Path(audit_mode) / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.%d%s' % (rank, step, suffix))
                    old = output / 'baseline' / relative
                    new = output / 'current' / relative
                    if integration.physical_payload(old) != integration.physical_payload(new):
                        raise ValueError('Default-off regression changed physical payload: ' + str(relative))
                    compared += 1
    result = dict(status='PASS', baseline_commit=commit, baseline_source='git archive HEAD lib bin tests/cbf',
                  temporary_helper_only_change='Replace git rev-parse HEAD expression with archived commit literal',
                  archived_source_tar_sha256=digest(archive), current_binary_sha256=original_hash,
                  baseline_binary_sha256=digest(baseline_build / 'CitcomSFull'), ranks=ranks,
                  volume_mode='default off in all four real MPI runs',
                  across_version_identical_physical_payloads=compared, independent_audit_checks=tests,
                  exclusions='Only unused native array index-zero slots; state sidecars compared exactly')
    (output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
