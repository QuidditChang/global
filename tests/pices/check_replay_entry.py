"""Compile original and adapted entry arithmetic with production input values.

No Python floating-point reimplementation of the C expressions: all three
functions below are extracted from actual C sources and use All_variables.
Usage: python3 check_replay_entry.py ORIGINAL_INPUT_CFG OUTPUT_DIR
"""
import argparse, json, os, subprocess
from pathlib import Path
from motion_window import parameters
from replay_entry import adapt, normalization


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    instructions = (root/'lib/Instructions.c').read_text()
    properties = (root/'module/setProperties.c').read_text()
    adapted = adapt(instructions, properties)
    # Fail closed if the C temporary no longer denotes the input radius.
    try:
        adapt(instructions.replace('E->data.radius_km = tmp / 1e3;',
                                   'E->data.radius_km = tmp / 1e3; tmp = 1;'), properties)
    except ValueError:
        pass
    else:
        raise AssertionError('clobbered radius accepted')
    functions = []
    for name, text in [('producer', properties), ('old_entry', instructions), ('fixed_entry', adapted)]:
        functions.append('void %s(struct All_variables *E,float radius) { float tmp=radius;\n%s\n}\n'
                         % (name, normalization(text).group()))
    code = '#include <stdio.h>\n#include <stdlib.h>\n#include <string.h>\n#include "global_defs.h"\n'
    code += ''.join(functions)
    code += r'''
int main(int argc,char **argv) {
 struct All_variables *v[3];float radius;int i;FILE *f=fopen(argv[1],"r");
 for(i=0;i<3;i++)v[i]=calloc(1,sizeof(*v[i]));
 if(!f || fscanf(f,"%f%f%f%f%f%f%f%f%f",&v[0]->data.rho0,&v[0]->data.g0,
  &v[0]->data.alpha0,&v[0]->data.ref_temperature,&radius,&v[0]->data.ref_viscosity,
  &v[0]->data.k0,&v[0]->data.Cp0,&v[0]->control.Atemp)!=9)return 2;
 fclose(f);
 v[0]->data.radius_km=radius/1e3;
 v[0]->data.kappa0=v[0]->data.k0/(v[0]->data.rho0*v[0]->data.Cp0);
 for(i=1;i<3;i++)memcpy(v[i],v[0],sizeof(*v[0]));
 producer(v[0],radius);old_entry(v[1],radius);fixed_entry(v[2],radius);
 printf("{\"producer_Ra\":%.17g,\"old_Ra\":%.17g,\"fixed_Ra\":%.17g,"
        "\"producer_deltaT\":%.17g,\"old_deltaT\":%.17g,\"fixed_deltaT\":%.17g}\n",
        (double)v[0]->control.Atemp,(double)v[1]->control.Atemp,(double)v[2]->control.Atemp,
        (double)v[0]->data.ref_temperature,(double)v[1]->data.ref_temperature,(double)v[2]->data.ref_temperature);
 return memcmp(v[0],v[2],sizeof(*v[0]))!=0;
}
'''
    (out/'entry_check.c').write_text(code)
    subprocess.run([os.environ.get('MPICC','mpicc'),'-O2','-std=gnu99','-I'+str(root/'lib'),
                    str(out/'entry_check.c'),'-o',str(out/'entry_check')],check=True)
    cfg = parameters(a.input.read_text())
    values = [cfg[k] for k in ('rho0','g0','alpha0')]
    values += [str(float(cfg['Tbottom'])-float(cfg['Ttop'])),cfg['radius'],cfg['refvisc'],cfg['k0'],cfg['Cp0'],'690018560']
    results = {}
    for name, delta, radius in [('production',values[3],values[4]),('nonintegral_radius',values[3],'6371000.5'),('inverse_deltaT','-1',values[4])]:
        case = values[:];case[3]=delta;case[4]=radius
        path = out/(name+'.txt');path.write_text(' '.join(case)+'\n')
        r = subprocess.run([str(out/'entry_check'),str(path)],capture_output=True,text=True,check=True)
        results[name] = json.loads(r.stdout)
    assert results['production']['producer_Ra'] != results['production']['old_Ra'],results
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
    print('PASS: original entry discrepancy reproduced; adapted entry byte-identical to producer arithmetic in all cases')


if __name__=='__main__':
    main()
