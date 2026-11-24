import os, runpy, sys
scripts = os.environ['SCRIPTDIR']
sys.path.insert(0, scripts)
import pices_prepare_restart
source = '0608be89a18479a0d22374eb9b395be940c49d86'
if os.environ.get('PICES_RESTART_SOURCE_COMMIT') != source:
    raise RuntimeError('Unexpected restart source')
pices_prepare_restart.SOURCES = set(pices_prepare_restart.SOURCES) | {source}
runpy.run_path(os.path.join(scripts, 'PostProc.lsf.py'), run_name='__main__')
