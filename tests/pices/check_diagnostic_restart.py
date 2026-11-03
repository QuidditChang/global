"""Compile actual metadata gate: permit only audited commit-only transitions."""
from pathlib import Path
import subprocess, tempfile
source=(Path(__file__).resolve().parents[2]/'lib/Pices_checkpoint.c').read_text()
helper=source[source.index('static int coverage_restart_compatible('):source.index('/* Explicit, audited transition')]
harness=r'''
#include <stdlib.h>
#include <string.h>
#include <assert.h>
'''+helper+r'''
int main(void) {
 char actual[256],expected[256];
 const char *old="63b5641310e619a4e74ef1c92a3918110368022b";
 const char *fresh="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
 #define META(dst,commit,heat) snprintf(dst,sizeof(dst),"{\"solver_commit\":\"%s\",\"RH\":%d}",commit,heat)
 META(actual,old,20);META(expected,fresh,20);
 unsetenv("PICES_RESTART_SOURCE_COMMIT");assert(!coverage_restart_compatible(actual,expected));
 setenv("PICES_RESTART_SOURCE_COMMIT",old,1);assert(coverage_restart_compatible(actual,expected));
 META(expected,fresh,10);assert(!coverage_restart_compatible(actual,expected));
 META(expected,fresh,20);actual[60]='X';assert(!coverage_restart_compatible(actual,expected));
 META(actual,fresh,20);assert(!coverage_restart_compatible(actual,expected));
 return 0;
}
'''
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'test.c').write_text('#include <stdio.h>\n'+harness)
 subprocess.run(['cc','-std=gnu11',str(p/'test.c'),'-o',str(p/'test')],check=True)
 subprocess.run([str(p/'test')],check=True)
print('PASS diagnostic restart metadata gate')
