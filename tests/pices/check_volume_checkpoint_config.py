"""Compile the real volume parameter reader and checkpoint implementation.

The checkpoint harness stubs only single-rank MPI collectives, not file formats,
hashes, metadata, state serialization or compatibility checks. Real 12/24-rank
restart tests additionally exercise the complete particle payload and solver.
"""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
COMPILER = os.environ.get('MPICC', '/usr/local/bin/mpicc')
COMMIT = '0123456789abcdef0123456789abcdef01234567'

PARAMETER_MAIN = r'''
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
#include "parsing.h"
void pices_fail(struct All_variables *E,const char *why) {
 fprintf(stderr,"PICES_ERROR: %s\n",why);exit(72);
}
void parallel_process_termination(void) {exit(73);}
void unique_copy_file(struct All_variables *E,char *name,char *comment) {}
PARAMETERS
int main(int argc,char **argv) {
 struct All_variables *E=calloc(1,sizeof(*E));E->fp=stdout;
 setup_parser(E,argv[1]);pices_parameters(E);
 printf("{\"mode\":%d,\"slot\":%d,\"hops\":%d,\"iterations\":%d,\"tolerance\":%.17g}\n",
 E->pices.volume_mode,E->pices.volume_slot,E->pices.volume_hops,
 E->pices.volume_iterations,E->pices.volume_tolerance);
 shutdown_parser(E);free(E);return 0;
}
'''

CHECKPOINT_MAIN = r'''
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
static int one_barrier(MPI_Comm world){return 0;}
static int one_allgather(const void *send,int n,MPI_Datatype a,void *recv,int m,MPI_Datatype b,MPI_Comm world) {
 assert(n==65 && m==65);memcpy(recv,send,65);return 0;
}
#define MPI_Barrier one_barrier
#define MPI_Allgather one_allgather
void pices_fail(struct All_variables *E,const char *why) {
 fprintf(stderr,"PICES_ERROR: %s\n",why);exit(72);
}
int pices_motion_enabled(struct All_variables *E){return 0;}
float find_age_in_MY(struct All_variables *E){return 0;}
#include "Pices_checkpoint.c"
static struct All_variables *model(int mode) {
 struct All_variables *E=calloc(1,sizeof(*E));int d;
 E->fp=stdout;E->parallel.nproc=E->parallel.nprocx=E->parallel.nprocy=E->parallel.nprocz=1;
 E->lmesh.nox=E->lmesh.noy=E->lmesh.noz=2;E->lmesh.nno=8;E->lmesh.nel=E->lmesh.npno=1;E->lmesh.neq=24;
 E->sphere.caps_per_proc=E->sphere.capid[1]=1;
 E->trace.nflavors=25;E->trace.number_of_basic_quantities=12;
 E->trace.number_of_extra_quantities=2+(mode!=0);E->trace.ntracers[1]=2;
 E->pices.enabled=E->pices.eba=1;E->pices.slot=1;E->pices.volume_mode=mode;
 E->pices.volume_slot=mode?2:-1;E->pices.volume_hops=2;E->pices.volume_iterations=1000;
 E->pices.volume_tolerance=1e-12;E->pices.volume_initialized=mode!=0;
 E->pices.volume_baseline[0]=.96;E->pices.volume_baseline[24]=.04;
 E->pices.volume_source[0]=-.01;E->pices.volume_source[1]=.01;
 E->composition.on=1;E->composition.ncomp=24;E->composition.ibuoy_type=1;
 E->composition.buoyancy_ratio=calloc(24,sizeof(double));
 E->control.RH=1;E->viscosity.RHEOL=1;E->viscosity.num_mat=1;E->advection.fixed_timestep=1;
 E->refstate.rho=calloc(3,sizeof(double));E->refstate.thermal_expansivity=calloc(3,sizeof(double));
 E->refstate.gravity=calloc(3,sizeof(double));E->refstate.heat_capacity=calloc(3,sizeof(double));
 for(d=1;d<=3;d++) {
  E->x[1][d]=calloc(9,sizeof(double));E->sphere.cap[1].TB[d]=calloc(9,sizeof(float));
  E->sphere.cap[1].V[d]=calloc(9,sizeof(float));E->sphere.cap[1].V[d][1]=d*.25f;
 }
 return E;
}
static void payload(struct All_variables *E,const char *path) {
 int h[8]={2,2,2,1,1,1,1,0},tr[5]={12,E->trace.number_of_extra_quantities,25,2,2};
 int cs[5]={0,0,0,0,24};float times[3]={0,0,0};
 long offset=44+16+16L*9+16+8+8L*(2+24);
 long bytes=offset+16+20+(8L*(6+tr[1])+4)*3+20+8L*24*4;
 unsigned char *data=calloc(bytes,1);FILE *f=fopen(path,"wb");
 memcpy(data,h,sizeof(h));memcpy(data+32,times,sizeof(times));
 memcpy(data+offset+16,tr,sizeof(tr));memcpy(data+bytes-20-8L*24*4,cs,sizeof(cs));
 assert(fwrite(data,1,bytes,f)==(size_t)bytes);assert(!fclose(f));free(data);
}
int main(int argc,char **argv) {
 const char *action=argv[1];int mode=atoi(argv[2]),f;
 struct All_variables *E=model(mode);char before[65],after[65];
 if(!strcmp(action,"publish") || !strcmp(action,"bad_sources")) {
  if(!strcmp(action,"bad_sources"))E->pices.volume_source[24]=.001;
  payload(E,"checkpoint.tmp");pices_checkpoint_publish(E,"checkpoint.tmp","checkpoint");return 0;
 }
 if(!strcmp(action,"preflight")){pices_checkpoint_preflight(E,"checkpoint");return 0;}
 if(!strcmp(action,"restore") || !strcmp(action,"changed_setting") || !strcmp(action,"changed_mode")) {
  memset(E->pices.volume_baseline,0,sizeof(E->pices.volume_baseline));
  memset(E->pices.volume_source,0,sizeof(E->pices.volume_source));E->pices.volume_initialized=0;
  pices_checkpoint_preflight(E,"checkpoint");pices_checkpoint_restore_state(E,"checkpoint");
  if(!strcmp(action,"changed_setting"))E->pices.volume_hops++;
  if(!strcmp(action,"changed_mode"))E->pices.volume_mode=3-mode;
  pices_checkpoint_check_state(E,"checkpoint");
  if(mode) {
   assert(E->pices.volume_initialized && E->pices.volume_baseline[0]==.96 && E->pices.volume_baseline[24]==.04);
   assert(E->pices.volume_source[0]==-.01 && E->pices.volume_source[1]==.01);
   for(f=2;f<24;f++)assert(E->pices.volume_source[f]==0 && E->pices.volume_baseline[f]==0);
  }
  return 0;
 }
 if(!strcmp(action,"fingerprint")) {
  fingerprint(E,before,0);E->pices.volume_hops++;fingerprint(E,after,0);
  assert((strcmp(before,after)!=0)==(mode!=0));return 0;
 }
 return 99;
}
'''


class VolumeCheckpointConfigTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='pices-volume-schema-')
        cls.build = Path(cls.temporary.name)
        source = (ROOT / 'lib/Pices.c').read_text()
        parameters = source[source.index('void pices_parameters('):source.index('void pices_validate(')]
        (cls.build / 'parameters.c').write_text(PARAMETER_MAIN.replace('PARAMETERS', parameters))
        (cls.build / 'checkpoint.c').write_text(CHECKPOINT_MAIN)
        flags = [COMPILER, '-std=gnu99', '-w', '-Wno-error=implicit-function-declaration',
                 '-Wno-error=implicit-int', '-I' + str(ROOT / 'lib'),
                 '-DPICES_SOLVER_COMMIT="' + COMMIT + '"']
        subprocess.run(flags + [str(cls.build / 'parameters.c'), str(ROOT / 'lib/Parsing.c'),
                                '-lm', '-o', str(cls.build / 'parameters')], check=True)
        subprocess.run(flags + [str(cls.build / 'checkpoint.c'), '-lm',
                                '-o', str(cls.build / 'checkpoint')], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def config(self, extra, expect=0):
        with tempfile.TemporaryDirectory() as folder:
            cfg = Path(folder) / 'case.cfg'
            cfg.write_text('energy_solver=pices\npices_eba=on\n' + extra)
            run = subprocess.run([str(self.build / 'parameters'), str(cfg)], capture_output=True, text=True)
            self.assertEqual(run.returncode, expect, run.stderr)
            return json.loads(run.stdout.splitlines()[-1]) if not expect else run.stderr

    def checkpoint(self, action, mode, directory, expect=0):
        run = subprocess.run([str(self.build / 'checkpoint'), action, str(mode)], cwd=directory,
                             capture_output=True, text=True)
        self.assertEqual(run.returncode, expect, run.stderr)
        return run

    def test_actual_standalone_parser_defaults_and_modes(self):
        self.assertEqual(self.config(''), dict(mode=0, slot=-1, hops=2, iterations=1000, tolerance=1e-12))
        for label, mode in (('off', 0), ('shadow', 1), ('active', 2)):
            self.assertEqual(self.config('pices_volume_mode=' + label + '\n')['mode'], mode)
        row = self.config('pices_volume_mode=active\npices_volume_hops=3\npices_volume_iterations=432\npices_volume_tolerance=5e-13\n')
        self.assertEqual((row['hops'], row['iterations'], row['tolerance']), (3, 432, 5e-13))

    def test_invalid_parameters_are_rejected_without_silent_clamping(self):
        for extra in ('pices_volume_mode=other', 'pices_volume_hops=0', 'pices_volume_hops=9',
                      'pices_volume_iterations=0', 'pices_volume_iterations=100001',
                      'pices_volume_hops=2.5', 'pices_volume_hops=words',
                      'pices_volume_iterations=1000abc', 'pices_volume_iterations=1e3',
                      'pices_volume_tolerance=nan', 'pices_volume_tolerance=-1',
                      'pices_volume_tolerance=1e-12abc', 'pices_volume_tolerance=words',
                      'pices_volume_tolerance=1e-7',
                      'pices_volume_mode=active\npices_eba=off',
                      'pices_volume_mode=active\nenergy_solver=pg'):
            with self.subTest(extra=extra):
                self.config(extra + '\n', 72)

    def test_modes_and_tunings_are_registered_in_pyre_and_C_bridge(self):
        inventory = (ROOT / 'CitcomS/Components/Advection_diffusion/Advection_diffusion.py').read_text()
        bridge = (ROOT / 'module/setProperties.c').read_text()
        for key, kind in (('mode', 'String'), ('hops', 'Int'), ('iterations', 'Int'), ('tolerance', 'Double')):
            name = 'pices_volume_' + key
            self.assertIn('"' + name + '"', inventory)
            self.assertIn('get' + kind + 'Property(properties, "' + name + '"', bridge)

    def test_authenticated_state_roundtrip_and_metadata(self):
        for mode in (0, 1, 2):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as folder:
                directory = Path(folder)
                self.checkpoint('publish', mode, directory)
                meta = json.loads((directory / 'checkpoint.pices.json').read_text())
                self.assertEqual(meta['schema'], 7 if mode else 3)
                self.assertEqual(len(meta['extraq']), 3 if mode else 2)
                self.assertEqual((directory / 'checkpoint.pices.state').stat().st_size,
                                 3 * 8 * 4 + 2 * 4 + (4 * 4 + 50 * 8 if mode else 0))
                if mode:
                    self.assertEqual(meta['extraq'][-1], dict(name='volume_weight', slot=2))
                    self.assertEqual(meta['composition_volume']['mode'], 'active' if mode == 2 else 'shadow')
                else:
                    self.assertNotIn('composition_volume', meta)
                self.checkpoint('restore', mode, directory)
                self.checkpoint('fingerprint', mode, directory)

    def test_old_weightless_checkpoint_cannot_start_new_mode(self):
        with tempfile.TemporaryDirectory() as folder:
            self.checkpoint('publish', 0, folder)
            for mode in (1, 2):
                result = self.checkpoint('preflight', mode, folder, 72)
                self.assertIn('original particle weights', result.stderr)

    def test_parameter_mode_and_payload_changes_cannot_restart_silently(self):
        with tempfile.TemporaryDirectory() as folder:
            self.checkpoint('publish', 2, folder)
            for action in ('changed_setting', 'changed_mode'):
                result = self.checkpoint(action, 2, folder, 72)
                self.assertIn('metadata/physics', result.stderr)
            self.checkpoint('preflight', 0, folder, 72)
            state = Path(folder) / 'checkpoint.pices.state'
            content = bytearray(state.read_bytes());content[-1] ^= 1;state.write_bytes(content)
            result = self.checkpoint('preflight', 2, folder, 72)
            self.assertIn('checksum mismatch', result.stderr)

    def test_invalid_conversion_ledger_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as folder:
            result = self.checkpoint('bad_sources', 2, folder, 72)
            self.assertIn('checkpoint sources', result.stderr)
            self.assertFalse((Path(folder) / 'checkpoint.pices.manifest').exists())


if __name__ == '__main__':
    unittest.main()
