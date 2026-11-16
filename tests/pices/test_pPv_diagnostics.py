"""Execute production pPv profile, HDF metadata and checkpoint hash code.

The fixtures supply geometry, horizontal averaging and HDF attribute sinks.
Phase fractions, profile getters, metadata selection and SHA serialization are
the actual production implementations. No Stokes phase cache is initialized.
"""
import importlib.util
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
COMPILER = os.environ.get('MPICC', '/usr/local/bin/mpicc')


def function(source, signature):
    start = source.index(signature)
    opening = source.index('{', start)
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[start:end]


PROFILE_FIXTURE = r'''
void parallel_process_termination(void) {abort();}
void remove_horiz_ave2(struct All_variables *E,double **f) {
 int z,i;for(z=1;z<=2;z++) {
  double avg=0;for(i=z;i<=8;i+=2)avg+=f[1][i]/4;
  for(i=z;i<=8;i+=2)f[1][i]-=avg;
 }
}
/* Existing k(T), k(depth), k(composition) are a fixed baseline here. */
double conductivity_temperature_factor(struct All_variables *E,double t){return 2;}
double conductivity_element_prefactor(struct All_variables *E,int m,int e,double k){return 3;}
'''

PROFILE_MAIN = r'''
int main(int argc,char **argv) {
 struct All_variables *E=calloc(1,sizeof(*E));int p,i,j,on=atoi(argv[1]);
 E->control.pPv_enabled=on;E->control.pPv_k_ratio=2;
 E->mesh.nsd=3;E->lmesh.noz=2;E->lmesh.nno=8;E->sphere.caps_per_proc=1;E->sphere.ro=1;
 E->control.Atemp=4;E->data.rho0=4000;E->data.ref_temperature=1000;
 E->data.alpha0=.001;E->data.k0=5;
 E->refstate.rho=calloc(3,sizeof(double));E->refstate.gravity=calloc(3,sizeof(double));
 E->refstate.Tref=calloc(3,sizeof(double));E->refstate.thermal_expansivity=calloc(3,sizeof(double));
 for(i=1;i<=2;i++){E->refstate.rho[i]=1;E->refstate.gravity[i]=1;
  E->refstate.Tref[i]=.5;E->refstate.thermal_expansivity[i]=1;}
 E->T[1]=calloc(9,sizeof(double));E->sx[1][3]=calloc(9,sizeof(double));
 E->ien[1]=calloc(2,sizeof(struct IEN));
 for(j=1;j<=3;j++)E->x[1][j]=calloc(9,sizeof(double));
 for(i=1;i<=8;i++){E->T[1][i]=.125*i;E->sx[1][3][i]=.5;E->x[1][3][i]=.5;
  E->ien[1][1].node[i]=i;E->N.vpt[GNVINDEX(i,1)]=.125;}
 for(p=0;p<PHASE_TRANSITIONS;p++) {
  E->control.phase[p].depth=p==pPv_PHASE_INDEX?.25:.125*(p+1);
  E->control.phase[p].transT=.5;E->control.phase[p].clapeyron=.25;
  E->control.phase[p].inv_width=8;E->control.phase[p].density_jump=40*(p+1);
  E->control.phase[p].Ra=p+1;
 }
 E->output.profile_flags[PROFILE_FLAG_BUOYANCY_TOTAL]=1;
 /* NULL phase_B pointers deliberately reproduce pre-Stokes restart output. */
 prepare_buoyancy_profile_fields(E);
 printf("%d %d %d %d\n",PROFILE_FLAG_BUOYANCY_TOTAL,
  profile_flag(E,PROFILE_FLAG_pPv_FRACTION),profile_flag(E,PROFILE_FLAG_BUOYANCY_PHASE_pPv),
  profile_flag(E,PROFILE_FLAG_BUOYANCY_TOTAL));
 printf("%.17g\n",conductivity_at_gp(E,1,1,1));
 for(i=1;i<=8;i++)printf("%.17g %.17g %.17g %.17g\n",
  density_total_at_node(E,1,i),pPv_fraction_profile_at_node(E,1,i),
  buoyancy_phase_cache[pPv_PHASE_INDEX][1][i],buoyancy_total_cache[1][i]);
 free_derived_profile_fields(E);return 0;
}
'''

HDF_FIXTURE = r'''
typedef int hid_t;
typedef int herr_t;
static int set_attribute_float_vector(int id,const char *key,int count,const float *v) {
 int i;printf("%s",key);for(i=0;i<count;i++)printf(" %.17g",(double)v[i]);puts("");return 0;
}
static int set_attribute_float(int id,const char *key,float v) {printf("%s %.17g\n",key,(double)v);return 0;}
static int set_attribute_int(int id,const char *key,int v) {printf("%s %d\n",key,v);return 0;}
'''

HASH_MAIN = r'''
int main(void) {
 struct All_variables *E=model(0);char h[65];int i;
 struct Phase_transition *p=&E->control.phase[pPv_PHASE_INDEX];
 float *values[]={&p->depth,&p->density_jump,&p->entropy_jump,&p->Ra,
  &p->clapeyron,&p->transT,&p->inv_width,&E->control.pPv_k_ratio,&E->control.pPv_eta_ratio};
 fingerprint(E,h,0);printf("legacy %s\n",h);
 fingerprint(E,h,1);printf("legacy_coverage %s\n",h);
 for(i=0;i<9;i++)*values[i]=(i+1)*.125f;
 E->control.pPv_latent_enabled=1;
 fingerprint(E,h,0);printf("disabled_populated %s\n",h);
 E->control.pPv_enabled=1;fingerprint(E,h,0);printf("enabled %s\n",h);
 for(i=0;i<9;i++) {
  float saved=*values[i];*values[i]+=.125f;
  fingerprint(E,h,0);printf("parameter_%d %s\n",i,h);*values[i]=saved;
 }
 E->control.pPv_latent_enabled=0;fingerprint(E,h,0);printf("latent_off %s\n",h);
 E->control.pPv_latent_enabled=1;E->control.pPv_enabled=0;
 E->control.phase[0].depth=.125f;fingerprint(E,h,0);printf("legacy_changed %s\n",h);
 return 0;
}
'''


class pPvDiagnostics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='pPv-diagnostics-')
        cls.build = Path(cls.temporary.name)
        profile = (ROOT/'lib/Profile_output.c').read_text()
        phase = (ROOT/'lib/Phase_change.c').read_text()
        includes = ('#include <math.h>\n#include <stdio.h>\n#include <stdlib.h>\n'
                    '#include "element_definitions.h"\n#include "global_defs.h"\n'
                    '#include "phase_change.h"\n')
        source = [includes, PROFILE_FIXTURE,
                  function(profile, 'enum Profile_flag {')+';',
                  profile[profile.index('static double **ala_residual_log10_abs_cache;'):
                          profile.index('extern void parallel_process_termination();')]]
        for signature in ('int phase_change_count(', 'void phase_change_state(',
                          'float phase_change_reference_fraction(',
                          'float phase_change_fraction_at_temperature(',
                          'double phase_change_density_fraction(',
                          'double pPv_fraction_at_state(', 'double pPv_fraction_at_node(',
                          'double pPv_fraction_at_gp(', 'double pPv_conductivity_factor('):
            source.append(function(phase, signature))
        for signature in ('static double density_total_at_node(',
                          'static double conductivity_at_gp(',
                          'static double pPv_fraction_profile_at_node(',
                          'static int profile_flag(', 'static void *profile_allocate(',
                          'static double **allocate_profile_field(',
                          'static void free_profile_field(', 'static int any_profile_flag(',
                          'static void prepare_buoyancy_profile_fields(',
                          'static void free_derived_profile_fields('):
            source.append(function(profile, signature))
        source.append(PROFILE_MAIN)
        cls.compile('profile', '\n'.join(source))

        hdf = function((ROOT/'lib/Output_h5.c').read_text(),
                       'void h5output_meta(struct All_variables *E)')
        declarations = '\n'.join(re.findall(r'    float phase_.*;', hdf))
        block = hdf[hdf.index('    for(phase_index=0;'):
                    hdf.index('    /*\n     * Solver.inventory')]
        hdf_main = r'''
int main(int argc,char **argv) {
 struct All_variables state={0},*E=&state;int i;E->control.pPv_enabled=atoi(argv[1]);
 E->control.pPv_latent_enabled=1;E->control.pPv_k_ratio=2;E->control.pPv_eta_ratio=.25;
 for(i=0;i<PHASE_TRANSITIONS;i++) {
  E->control.phase[i].depth=.125*(i+1);E->control.phase[i].density_jump=40*(i+1);
  E->control.phase[i].entropy_jump=-.125;E->control.phase[i].inv_width=8;
 }
 metadata(E);return 0;
}
'''
        cls.compile('hdf', includes+HDF_FIXTURE+
                    'static void metadata(struct All_variables *E) {\n'
                    'int phase_index,input=0,status;\n'+declarations+'\n'+block+'}\n'+hdf_main)

        spec = importlib.util.spec_from_file_location('checkpoint_fixture',
                            ROOT/'tests/pices/check_volume_checkpoint_config.py')
        fixture = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fixture)
        checkpoint = fixture.CHECKPOINT_MAIN.split('int main(int argc,char **argv)')[0]
        cls.compile('hash', checkpoint+HASH_MAIN)

    @classmethod
    def compile(cls, name, source):
        path = cls.build/(name+'.c')
        path.write_text(source)
        subprocess.run([COMPILER, '-std=gnu99', '-w', '-I'+str(ROOT/'lib'),
                        str(path), '-lm', '-o', str(cls.build/name)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def run_fixture(self, name, *args):
        return subprocess.check_output([str(self.build/name), *map(str,args)], text=True)

    def profile(self, on):
        rows = self.run_fixture('profile', on).splitlines()
        return list(map(int,rows[0].split())), float(rows[1]), [list(map(float,r.split())) for r in rows[2:]]

    def test_profiles_reconstruct_density_and_buoyancy_before_stokes(self):
        off_flags,off_k,off = self.profile(0)
        on_flags,on_k,on = self.profile(1)
        self.assertEqual(off_flags,[22,0,0,1])
        self.assertEqual(on_flags,[22,1,1,1])
        fractions=[.5*(1+math.tanh(8*(.25-.25*(.125*i-.5)))) for i in range(1,9)]
        for i,(a,b,x) in enumerate(zip(off,on,fractions)):
            self.assertEqual(a[1],0)
            self.assertAlmostEqual(b[1],x,places=14)
            self.assertAlmostEqual(b[0]-a[0],160*x,places=4)  # nodal X is stored as float
            mean=sum(fractions[i%2::2])/4
            expected=-4*(x-mean)
            self.assertAlmostEqual(b[2],expected,places=6)
            self.assertAlmostEqual(b[3]-a[3],b[2],places=12)

    def test_conductivity_uses_one_pPv_factor_at_gp_state(self):
        _,off,_=self.profile(0)
        _,on,_=self.profile(1)
        x=.5*(1+math.tanh(8*(.25-.25*(.5625-.5))))
        self.assertEqual(off,30)
        self.assertAlmostEqual(on,30*(1+x),places=12)

    def test_hdf_legacy_vectors_stay_three_entries_and_pPv_is_separate(self):
        def attributes(on):
            return {row.split()[0]:list(map(float,row.split()[1:]))
                    for row in self.run_fixture('hdf',on).splitlines()}
        off,on=attributes(0),attributes(1)
        self.assertEqual(len(off),6)
        self.assertTrue(all(len(v)==3 for v in off.values()))
        for key in off:
            self.assertEqual(off[key],on[key])
        self.assertEqual(set(on)-set(off),{
            'pPv_enabled','pPv_latent_enabled','pPv_depth','pPv_delta_rho',
            'pPv_delta_s','pPv_Ra','pPv_width','pPv_clapeyron','pPv_transT',
            'pPv_k_ratio','pPv_eta_ratio'})
        self.assertEqual(on['pPv_depth'],[.5])
        self.assertEqual(on['pPv_width'],[.125])

    def test_disabled_fingerprint_matches_frozen_af0605d_fixture(self):
        if sys.byteorder!='little':
            self.skipTest('Frozen native-byte baseline was generated on a little-endian host')
        hashes=dict(row.split() for row in self.run_fixture('hash').splitlines())
        # Baselines executed from af0605d with the same model(0) fixture and
        # legacy PHASE_TRANSITIONS=3; both coverage-policy paths are retained.
        self.assertEqual(hashes['legacy'],'1199283228beca19c7dd26b22a03cf40c5ac23882855c861397e164e7ef6fb1a')
        self.assertEqual(hashes['legacy_coverage'],'b704a56c3af0846ae9a2756c93e0b27118af4fe64c7bc2910a368dcf85990c13')
        self.assertEqual(hashes['disabled_populated'],hashes['legacy'])
        self.assertNotEqual(hashes['legacy_changed'],hashes['legacy'])

    def test_enabled_fingerprint_covers_every_pPv_setting(self):
        hashes=dict(row.split() for row in self.run_fixture('hash').splitlines())
        self.assertNotEqual(hashes['enabled'],hashes['legacy'])
        for key in ['latent_off']+['parameter_'+str(i) for i in range(9)]:
            self.assertNotEqual(hashes[key],hashes['enabled'],key)


if __name__=='__main__':
    unittest.main()
