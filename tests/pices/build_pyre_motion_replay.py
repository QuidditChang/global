"""Relink the existing embedded Pyre build into an isolated replay executable.

Reuse the producer's compiled initializer, mover and Python bindings. Replace
only a COPY of its checkpoint archive member, adding diagnostics and a replay
exit after strict checkpoint acceptance. No make/install in the producer tree.
"""
import argparse, hashlib, json, os, re, shlex, shutil, subprocess
from pathlib import Path
from motion_window import SOLVER_COMMIT
from replay_checkpoint_instrument import instrument, VERSION

ROOT = Path(__file__).resolve().parents[2]
COMPILE = '$(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(libCitcomS_a_CFLAGS) $(CFLAGS)'
LINK = ('$(PYLINK) $(PYTHON_LDFLAGS) $(PYTHON_LINKFORSHARED) '
        '$(mpipycitcoms_LDFLAGS) $(mpipycitcoms_OBJECTS) $(mpipycitcoms_LDADD) '
        '$(PYTHON_EGG_LDFLAGS) $(PYTHON_BLDLIBRARY) $(PYTHON_LIBS) '
        '$(PYTHON_MODLIBS) $(PYTHON_SYSLIBS) $(LIBS) $(PYTHON_LDLAST)')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()


def run_logged(args, cwd, path):
    with path.open('w') as log:
        r=subprocess.run(args,cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError('Build failed; '+str(path)+'\n'+'\n'.join(path.read_text(errors='replace').splitlines()[-50:]))


def make_values(directory, output, expressions):
    probe = output/(directory.name+'_variables.mk')
    probe.write_text('.PHONY: pices_read_variables\npices_read_variables:\n' +
                     ''.join('\t$(info PICES_VAR_%s=%s)\n'%(k,v) for k,v in expressions.items())+'\t@:\n')
    # -n expands recipes but never compiles, relinks or installs anything.
    r = subprocess.run(['make','--no-print-directory','-n','-f','Makefile','-f',str(probe),
                        'pices_read_variables'],cwd=directory,capture_output=True,text=True,check=True)
    values = dict(line[len('PICES_VAR_'):].split('=',1) for line in r.stdout.splitlines() if line.startswith('PICES_VAR_'))
    if values.keys()!=expressions.keys(): raise ValueError('Cannot read producer make variables: '+r.stdout+r.stderr)
    return values


def checkpoint_source(source):
    text = instrument(source)
    text = text.replace('#include "replay_checkpoint_diagnostic.h"',
                        '#include "replay_checkpoint_diagnostic.h"\n#include "motion_replay_driver.h"')
    start = text.index('void pices_checkpoint_check_state(')
    end = text.index('/* U is the Stokes iterate;',start)
    part = text[start:end]
    closing = part.rfind('\n}')
    if closing<0: raise ValueError('Checkpoint acceptance boundary changed')
    part = part[:closing] + '''
 if(E->parallel.me==0)fprintf(stderr,"PICES_REPLAY_ENTRY entry=producer_objects checkpoint=accepted\\n");
 motion_replay(E);
 /* All ranks stop before the Python Controller can save or advance physics. */
 MPI_Barrier(E->parallel.world);
 MPI_Finalize();
 exit(0);
''' + part[closing:]
    return text[:start]+part+text[end:]


def build(producer, output, python_home=None):
    producer=producer.resolve();output=output.resolve()
    if output.exists(): raise ValueError('Use a new isolated build directory: '+str(output))
    stamp=(producer/'pices-p0-build.commit').read_text().strip()
    paths=['lib','bin','module','CitcomS','config_script','configure.ac','Makefile.am','m4',
           ':(exclude,glob)**/Makefile.in',':(exclude,glob)**/*.pyc']
    for args in ([SOLVER_COMMIT,'--',*paths],[SOLVER_COMMIT,stamp,'--',*paths]):
        subprocess.run(['git','-C',str(ROOT),'diff','--exit-code',*args],check=True)
    archive=producer/'lib/libCitcomS.a';module=producer/'module/libCitcomSLibmodule.a'
    binary=producer/'install/bin/mpipycitcoms'
    site=producer/'install/lib/python2.6/site-packages/CitcomS'
    for p in (archive,module,binary,producer/'lib/Makefile',producer/'bin/Makefile',site/'Solver/Solver.py'):
        if not p.is_file(): raise ValueError('Missing embedded production build artifact: '+str(p))
    original={str(p):sha(p) for p in (archive,module,binary)}
    output.mkdir(parents=True)
    shutil.copytree(site,output/'PyreSite/CitcomS',ignore=shutil.ignore_patterns('*.pyc','*.pyo','__pycache__'))
    libvars=make_values(producer/'lib',output,dict(COMPILE=COMPILE,CC='$(CC)',AR='$(AR)',RANLIB='$(RANLIB)'))
    binvars=make_values(producer/'bin',output,dict(LINK=LINK,OBJECTS='$(mpipycitcoms_OBJECTS)',LDADD='$(mpipycitcoms_LDADD)'))
    ar=shlex.split(libvars['AR']);members=subprocess.check_output([*ar,'t',str(archive)],text=True).splitlines()
    matches=[s for s in members if re.fullmatch(r'(libCitcomS_a-)?Pices_checkpoint\.o',s)]
    if len(matches)!=1: raise ValueError('Expected one checkpoint object in producer archive')
    member=matches[0]
    (output/'Pices_checkpoint.c').write_text(checkpoint_source((ROOT/'lib/Pices_checkpoint.c').read_text()))
    for name in ('replay_checkpoint_diagnostic.h','motion_replay_driver.h'):
        shutil.copy2(ROOT/'tests/pices'/name,output/name)
    compile_args=shlex.split(libvars['COMPILE'])
    compile_args=[s for s in compile_args if not s.startswith('-DPICES_SOLVER_COMMIT=')]
    compile_args+=['-DPICES_SOLVER_COMMIT="'+SOLVER_COMMIT+'"','-c',str(output/'Pices_checkpoint.c'),'-o',str(output/member)]
    run_logged(compile_args,producer/'lib',output/'compile.log')
    copied=output/'libCitcomS.a';shutil.copy2(archive,copied)
    subprocess.run([*ar,'r',str(copied),member],cwd=output,check=True)
    subprocess.run([*shlex.split(libvars['RANLIB']),str(copied)],check=True)
    # Prove the archive splice changed exactly the intended member.
    if subprocess.check_output([*ar,'t',str(copied)],text=True).splitlines()!=members:
        raise ValueError('Archive member inventory changed')
    for name in members:
        # BSD ar lists its regenerated symbol index as a pseudo-member.
        if name==member or name.startswith('__.SYMDEF'): continue
        old=subprocess.check_output([*ar,'p',str(archive),name])
        new=subprocess.check_output([*ar,'p',str(copied),name])
        if old!=new: raise ValueError('Unexpected changed production object: '+name)
    link_args=shlex.split(binvars['LINK']);replacements=0
    for i,arg in enumerate(link_args):
        if arg=='-o': link_args[i+1]=str(output/'PICESMotionReplay')
        if (producer/'bin'/arg).resolve()==archive:
            link_args[i]=str(copied);replacements+=1
    if replacements!=1 or link_args.count('-o')!=1: raise ValueError('Unexpected producer link recipe')
    run_logged(link_args,producer/'bin',output/'link.log')
    with (output/'PICESMotionReplay').open('rb') as f:
        if f.read(2)==b'#!':raise ValueError('Link produced a libtool wrapper, not a portable embedded executable; see link.log')
    for p,expected in original.items():
        if sha(Path(p))!=expected: raise ValueError('Producer artifact changed during isolated build: '+p)
    objects={str((producer/'bin'/p).resolve()):sha((producer/'bin'/p).resolve())
             for p in shlex.split(binvars['OBJECTS']) if p.endswith('.o')}
    versions={}
    for flag in ('--version','-show'):
        r=subprocess.run([*shlex.split(libvars['CC']),flag],cwd=producer/'lib',capture_output=True,text=True)
        versions[flag]=dict(returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
    bootstrap=None
    if python_home is not None:
        python_home=python_home.resolve()
        paths=[output/'PyreSite',site.parent,python_home/'lib/python2.6/site-packages']
        for folder in (producer/'install/python',producer/'python',site.parent):paths+=sorted(folder.glob('*.egg'))
        env=dict(os.environ,PYTHONHOME=str(python_home),PYTHONPATH=':'.join(map(str,paths)),
                 PYTHONNOUSERSITE='1',I_MPI_FABRICS='shm')
        code='import CitcomSLib; from CitcomS.SimpleApp import SimpleApp; print("PICES_PYRE_BOOTSTRAP_OK")'
        r=subprocess.run([str(output/'PICESMotionReplay'),'-c',code],env=env,cwd=output,capture_output=True,text=True,timeout=90)
        (output/'pyre_bootstrap.log').write_text(r.stdout+r.stderr)
        if r.returncode or 'PICES_PYRE_BOOTSTRAP_OK' not in r.stdout:
            raise ValueError('Native Pyre bootstrap failed; '+str(output/'pyre_bootstrap.log')+'\n'+(r.stdout+r.stderr)[-4000:])
        bootstrap='PASS'
    manifest=dict(solver_commit=SOLVER_COMMIT,driver_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
                  isolated_compile=True,entry='producer_pyre_objects_v1',checkpoint_diagnostic_version=VERSION,
                  producer_build_commit=stamp,producer_artifacts=original,producer_entry_objects=objects,
                  replaced_archive_member=member,unchanged_archive_members=sum(n!=member and not n.startswith('__.SYMDEF') for n in members),
                  pyre_files={str(p.relative_to(output/'PyreSite')):sha(p) for p in (output/'PyreSite').rglob('*') if p.is_file()},
                  pyre_bootstrap=bootstrap,
                  compile_argv=compile_args,link_argv=link_args,compiler=versions,
                  checkpoint_source_sha256=sha(output/'Pices_checkpoint.c'),helper_sha256=sha(output/'motion_replay_driver.h'),
                  binary_sha256=sha(output/'PICESMotionReplay'))
    (output/'replay_build.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(output/'PICESMotionReplay')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path);p.add_argument('--producer-build',type=Path,default=ROOT)
    p.add_argument('--python-home',type=Path,required=True)
    a=p.parse_args();build(a.producer_build,a.output,a.python_home)
