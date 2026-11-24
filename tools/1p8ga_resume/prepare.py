"""Resume the same PICES run with its verified frozen runtime and unchanged physics."""
import sys
import os
sys.path.insert(0, os.environ.get("SCRIPTDIR", "/work/ess-zhangkd/CitcomS/scripts"))
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace

from _postproc_cleanup import Cleanup, digest
from pices_prepare_restart import config, rewrite, resume_online


def verify_runtime(root, runtime=None):
    runtime = root / 'Runtime' if runtime is None else Path(runtime)
    if root.resolve() not in runtime.resolve().parents:
        raise ValueError('Runtime directory outside run')
    manifest = json.loads((runtime / 'manifest.json').read_text())
    for record in manifest['copies']:
        destination = Path(record['destination'])
        if runtime.resolve() not in destination.resolve().parents:
            raise ValueError('Runtime destination outside run')
        for relative, sha in record['sha256'].items():
            p = destination if relative == '.' else destination / relative
            if digest(p) != sha:
                raise ValueError('Frozen runtime checksum mismatch: ' + str(p))
    if not manifest['copies']:
        raise ValueError('Empty runtime inventory')
    return manifest['solver_commit']


def prepare(root, template, storage, tag, runtime=None, compatible_source_commit=None):
    root = Path(os.path.abspath(str(root)))  # retain the rollover symlink
    if not root.is_dir():
        raise ValueError('Existing run directory required')
    if compatible_source_commit not in (None, '0608be89a18479a0d22374eb9b395be940c49d86'):
        raise ValueError('Unaudited checkpoint source commit')
    if (storage / 'request.json').exists():
        raise ValueError('Unfinished migration request: inspect it before resuming')
    old, wanted = config(root/'input.cfg'), config(template)
    if old.get('CitcomS.solver.tsolver', 'energy_solver') != 'pices':
        raise ValueError('Expected a PICES run')
    if old.getfloat('CitcomS.solver.tsolver', 'pices_cfl', fallback=0.24) != wanted.getfloat('CitcomS.solver.tsolver', 'pices_cfl', fallback=0.24):
        raise ValueError('This resume must preserve CFL')
    mesh = 'CitcomS.solver.mesher'
    ranks = old.getint(mesh,'nproc_surf',fallback=12)*old.getint(mesh,'nprocx')*old.getint(mesh,'nprocy')*old.getint(mesh,'nprocz')
    commit = verify_runtime(root, runtime)
    w = SimpleNamespace(root=root, ranks=ranks, prefix=old.get('CitcomS.solver','datafile'))
    cleaner = Cleanup(w, old)
    steps = [39400]  # Exact audited restart point; never silently roll back farther.
    selected = None
    for step in steps:
        if step <= 0: continue
        try:
            cleaner.validate_checkpoint(step)
            for rank in range(ranks):
                meta=json.loads((root/'DATA'/str(rank)/('%s.chkpt.%d.%d.pices.json'%(w.prefix,rank,step))).read_text())
                if meta['solver_commit'] not in (commit, compatible_source_commit) or meta.get('magic') != 'CITCOMS_EBA_PICES':
                    raise ValueError('Checkpoint does not match frozen runtime')
        except FileNotFoundError:
            print('SKIP incomplete checkpoint', step, flush=True)
            continue
        selected=step
        break
    if selected is None: raise ValueError('No complete positive checkpoint')
    if selected >= wanted.getint('CitcomS','steps'): raise ValueError('Already at step limit')
    for rank in range(ranks):
        if not (root/'DATA'/str(rank)/('%s.coord.%d'%(w.prefix,rank))).is_file():
            raise ValueError('Missing original coordinate rank=%d'%rank)
    backup=root/'PostProc/Online'/('resume-'+tag)
    backup.mkdir()  # refuse duplicate preparation
    for relative in ('input.cfg','restart_manifest.json','submitted.cfg','submitted.lsf','solver_commit.txt','runs_commit.txt','scripts_commit.txt','input.sha256','PostProc/Online/state.json','PostProc/Online/launch.json','PostProc/Online/solver.finished.json','PostProc/Online/result.json'):
        p=root/relative
        if p.exists():
            dst=backup/relative; dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(p,dst)
    text=rewrite((root/'input.cfg').read_text(),{
        ('CitcomS','steps'):wanted.get('CitcomS','steps'),
        ('CitcomS.job','walltime'):wanted.get('CitcomS.job','walltime'),
        ('CitcomS.solver','datadir'):str(root/'DATA/%RANK'),
        ('CitcomS.solver','datadir_old'):str(root/'DATA/%RANK'),
        ('CitcomS.solver.ic','restart'):'on',
        ('CitcomS.solver.ic','solution_cycles_init'):str(selected),
        ('CitcomS.solver.ic','zero_elapsed_time'):'off',
    })
    (root/'input.cfg').write_text(text)
    resume_online(root,selected,text)
    record=dict(step=selected,solver_commit=commit,logical_root=str(root),physical_root=str(root.resolve()),
                checkpoint_hash_verification='all_rank_payload_state_and_collective_manifest',
                compatible_source_commit=compatible_source_commit,physics='original_input_cfg',runtime=str(runtime or root/'Runtime'),backup=str(backup))
    (root/'restart_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)
    return record



SOURCE = '0608be89a18479a0d22374eb9b395be940c49d86'

def snapshot(root, code):
    import subprocess
    commit = subprocess.check_output(['git','-C',str(code),'rev-parse','HEAD'], universal_newlines=True).strip()
    if (code/'pices-p0-build.commit').read_text().strip() != commit:
        raise ValueError('Build stamp differs from source HEAD')
    subprocess.run(['git','-C',str(code),'diff','--exit-code','HEAD','--','lib','CitcomS','bin',
                    ':(exclude,glob)**/Makefile.in',':(exclude,glob)**/*.pyc'],check=True)
    if verify_runtime(root) != SOURCE:
        raise ValueError('Original frozen Runtime is not the audited producer')
    runtime = root / ('Runtime_host_fix_' + commit[:12])
    if runtime.exists():
        if verify_runtime(root, runtime) != commit:
            raise ValueError('Existing fixed runtime identity mismatch')
        return runtime
    runtime.mkdir()
    def inventory(folder):
        return {str(f.relative_to(folder)): digest(f) for f in sorted(folder.rglob('*'))
                if f.is_file() and f.suffix not in ('.pyc','.pyo') and not f.name.startswith('._')}
    src = code/'install'
    before = inventory(src)
    shutil.copytree(src, runtime/'install', ignore=shutil.ignore_patterns('*.pyc','*.pyo','__pycache__','._*'))
    if inventory(runtime/'install') != before or inventory(src) != before:
        raise ValueError('Build changed during snapshot')
    if (code/'pices-p0-build.commit').read_text().strip() != commit:
        raise ValueError('Build changed during snapshot')
    # The original dependency eggs stay in the original verified Runtime.
    (runtime/'eggs.path').write_text((root/'Runtime/eggs.path').read_text())
    (runtime/'manifest.json').write_text(json.dumps(dict(solver_commit=commit,copies=[
        dict(destination=str(runtime/'install'),sha256=before)]),indent=2)+'\n')
    return runtime

if __name__ == '__main__':
    import subprocess
    root=Path(os.environ['RUNROOT'])
    code=Path(os.environ['CODE_DIR'])
    storage=Path(os.environ['STORAGE_CONTROL'])
    tag=os.environ['LSB_JOBID']
    if not tag.isdigit(): raise ValueError('Expected LSF job ID')
    with (root/'PostProc/Online/watch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        runtime=snapshot(root,code)
        cfg=config(root/'input.cfg')
        if cfg.getint('CitcomS','steps') != 300000:
            raise ValueError('Unexpected step limit; expected original 300000')
        # Validate and preserve the restart payload before any mutation.
        w=SimpleNamespace(root=root,ranks=384,prefix='global')
        Cleanup(w,cfg).validate_checkpoint(39400)
        backup=root/'Restart'/'host_fix_39400'
        if backup.exists(): raise ValueError('Restart backup already exists; inspect previous attempt')
        for rank in range(384):
            meta=json.loads((root/'DATA'/str(rank)/('global.chkpt.%d.39400.pices.json'%rank)).read_text())
            if meta['solver_commit'] != SOURCE: raise ValueError('Unexpected checkpoint producer')
        for rank in range(384):
            target=backup/str(rank);target.mkdir(parents=True)
            for suffix in ('','.pices.json','.pices.state','.pices.manifest'):
                f=root/'DATA'/str(rank)/('global.chkpt.%d.39400%s'%(rank,suffix))
                shutil.copy2(f,target/f.name)
        prepare(root,root/'input.cfg',storage,tag,runtime,SOURCE)
        (root/'PostProc/Online'/('resume-'+tag)/'runtime.path').write_text(str(runtime)+'\n')
