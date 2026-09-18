"""Frozen launch isolation and fail-closed interpretation of HPC evidence."""
import contextlib
import ctypes
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT.parents[1]/"runs"
spec = importlib.util.spec_from_file_location("frozen", ROOT/"tools/strict_ala_frozen_current.py")
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)


def write_csv(path, records):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


class FrozenCurrentTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)/"experiment"
        with contextlib.redirect_stdout(io.StringIO()):
            frozen.prepare(RUNS, ROOT, self.root, False)

    def tearDown(self):
        self.temp.cleanup()

    def evidence(self):
        for case in frozen.CASES:
            folder = self.root/case
            data = folder/"DATA/0"
            (folder/"exit_status.txt").write_text("0\n")
            (data/"global_AhatP60_test.log").write_text(
                "STRICT_ALA_FROZEN_CURRENT phase=begin physical=abc viscosity=def matched=1 step=0 time=0\n"
                "STRICT_ALA_FROZEN_CURRENT phase=end physical=abc viscosity=def matched=1 step=0 time=0\n")
            if case == "ADJOINT":
                write_csv(data/"global.strict_ala_stage_B_adjoint.csv",
                          [dict(delta_scale=1e-13,repeatability_relative=0)]*42)
                write_csv(data/"global.strict_ala_stage_B_gauge.csv",
                          [dict(achieved_first=1e-9,achieved_second=1e-9,qT_S_gamma_q=1,
                                BT_split_relative=0,repeat_action_relative=0)]*4)
            else:
                write_csv(data/"global.strict_ala_stage_C_iterations.csv",
                          [dict(iteration=60,final_iterate=1,continuity_relative=.02,
                                momentum_relative=1e-5,krylov_drift=1)])
                write_csv(data/"global.strict_ala_stage_C_inner_solves.csv",
                          [dict(status="CONVERGED",seconds=1)])

    def analyze(self):
        with contextlib.redirect_stdout(io.StringIO()):
            frozen.analyze(self.root)
        return json.loads((self.root/"analysis.json").read_text())

    def test_controlled_differences_and_original_inputs_preserved(self):
        baseline = frozen.cfg_read(self.root/"BASE/case.cfg")
        for case in ("TIGHT","BPI","UNSCALED"):
            cfg = frozen.cfg_read(self.root/case/"case.cfg")
            differences = {k for k in cfg if cfg[k] != baseline.get(k)
                           and k[1] not in ("datadir","datadir_old")}
            self.assertEqual(differences,set(frozen.VARIANTS[case]))
        for case in frozen.CASES:
            for name in frozen.INPUTS[1:]:
                self.assertEqual(frozen.digest(RUNS/name),frozen.digest(self.root/case/name))
        self.assertEqual(baseline["CitcomS","steps"],"1")

    def test_refuses_existing_directory(self):
        with self.assertRaises(ValueError):
            frozen.prepare(RUNS, ROOT, self.root, False)

    def test_hpc_preflight_explains_missing_build_receipt(self):
        code=Path(self.temp.name)/"unstamped-code"
        code.mkdir()
        with self.assertRaisesRegex(ValueError, "stamp-build"):
            frozen.prepare(RUNS, code, Path(self.temp.name)/"hpc", True)

    def test_nonconvergence_is_valid_evidence_not_acceptance(self):
        self.evidence()
        result=self.analyze()
        self.assertTrue(result["valid"],result)
        self.assertFalse(result["cases"]["BASE"]["joint_converged"])

    def test_missing_guard_invalidates_comparison(self):
        self.evidence()
        (self.root/"BASE/DATA/0/global_AhatP60_test.log").write_text("")
        self.assertFalse(self.analyze()["valid"])

    def test_only_final_iterate_can_be_accepted(self):
        self.evidence()
        path=self.root/"BASE/DATA/0/global.strict_ala_stage_C_iterations.csv"
        records=frozen.rows(path)
        earlier=dict(records[0],iteration="59",final_iterate="0",continuity_relative="0.001")
        write_csv(path,[earlier]+records)
        result=self.analyze()
        self.assertTrue(result["valid"])
        self.assertFalse(result["cases"]["BASE"]["joint_converged"])

    def test_changed_physics_invalidates_comparison(self):
        self.evidence()
        path=self.root/"TIGHT/DATA/0/global_AhatP60_test.log"
        path.write_text(path.read_text().replace("physical=abc","physical=123"))
        self.assertFalse(self.analyze()["valid"])

    def test_mutated_input_invalidates_comparison(self):
        self.evidence()
        path=self.root/"BPI/refstate_ALA_strict.txt"
        path.write_text(path.read_text()+"# changed\n")
        self.assertFalse(self.analyze()["valid"])

    def test_missing_or_nonfinite_trajectory_invalidates_comparison(self):
        self.evidence()
        path=self.root/"BASE/DATA/0/global.strict_ala_stage_C_iterations.csv"
        path.write_text(path.read_text().replace("0.02","nan"))
        self.assertFalse(self.analyze()["valid"])

    def test_failed_adjoint_is_flagged_before_solver_attribution(self):
        self.evidence()
        path=self.root/"ADJOINT/DATA/0/global.strict_ala_stage_B_adjoint.csv"
        path.write_text(path.read_text().replace("1e-13","0.1"))
        self.assertIn("OPERATOR_OR_MPI_ADJOINT_FAILURE",self.analyze()["interpretation"])


class FrozenGuardKernelTest(unittest.TestCase):
    """Run the actual C state/hash guard with a serial reduction fixture."""
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        directory=Path(cls.temp.name)
        source=(ROOT/"lib/Drive_solvers.c").read_text()
        def extract(signature):
            start=source.index(signature)
            pos=source.index("{",start)+1
            depth=1
            while depth:
                depth+=(source[pos]=="{")-(source[pos]=="}")
                pos+=1
            return source[start:pos]
        prefix=r'''
#include <string.h>
#include <setjmp.h>
#include "element_definitions.h"
#include "global_defs.h"
static jmp_buf failed;
void myerror(struct All_variables *E,char *message) { longjmp(failed,1); }
static int serial_reduce(const void *a,void *b,int n,MPI_Datatype type,
                         MPI_Op op,MPI_Comm comm) {
    memcpy(b,a,n*(type==MPI_INT ? sizeof(int) : sizeof(unsigned long long)));
    return 0;
}
#define MPI_Allreduce serial_reduce
'''
        main=r'''
int run_guard(int mode) {
    static struct All_variables state;
    struct All_variables *E=&state;
    double t[2]={0,0.4},tdot[2]={0,0};
    float coarse[9]={0},fine[9]={0};
    int result;
    memset(E,0,sizeof(*E));
    E->sphere.caps_per_proc=1; E->lmesh.nno=1;
    E->mesh.nsd=3; E->mesh.gridmin=0; E->mesh.levmax=1;
    E->lmesh.NEL[0]=E->lmesh.NEL[1]=1;
    E->T[1]=t; E->Tdot[1]=tdot;
    E->EVI[0][1]=coarse; E->EVI[1][1]=fine;
    E->fp=fopen("/dev/null","w");
    setenv("STRICT_ALA_FROZEN_CURRENT","1",1);
    result=setjmp(failed);
    if(!result) {
        strict_ala_frozen_current_guard(E,1);
        if(mode==1) t[1]+=0.01;
        if(mode==2) coarse[8]=1;
        if(mode==3) fine[1]=1;
        if(mode==4) E->monitor.elapsed_time+=1;
        strict_ala_frozen_current_guard(E,0);
    }
    fclose(E->fp);
    unsetenv("STRICT_ALA_FROZEN_CURRENT");
    return result;
}
'''
        code=prefix+"\n"+"\n".join(extract(sig) for sig in (
            "static unsigned long long strict_ala_vc1_hash_bytes(",
            "static unsigned long long strict_ala_vc1_physical_state_hash(",
            "void strict_ala_frozen_current_guard("))+main
        (directory/"guard.c").write_text(code)
        subprocess.run([os.environ.get("MPICC","mpicc"),"-std=gnu99","-shared","-fPIC",
                        "-I"+str(ROOT/"lib"),str(directory/"guard.c"),
                        "-o",str(directory/"guard.so")],check=True,capture_output=True)
        cls.library=ctypes.CDLL(str(directory/"guard.so"))
        cls.library.run_guard.argtypes=[ctypes.c_int]
        cls.library.run_guard.restype=ctypes.c_int

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_unchanged_state_passes(self):
        self.assertEqual(self.library.run_guard(0),0)

    def test_temperature_clock_and_each_viscosity_level_are_protected(self):
        for mode in range(1,5):
            self.assertEqual(self.library.run_guard(mode),1)


if __name__ == "__main__":
    unittest.main()
