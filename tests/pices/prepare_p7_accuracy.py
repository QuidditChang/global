"""Matched-time, nested-grid follow-up to the unresolved P7 temperature gap."""
import argparse
import json
import re
from pathlib import Path


def edit(text, options):
    for key, value in options.items():
        pattern = '^' + re.escape(key) + '=.*$'
        line = key + '=' + str(value)
        text = re.sub(pattern, lambda _: line, text, flags=re.M) if re.search(pattern, text, re.M) else text + '\n' + line + '\n'
    return text


def prepare(runs):
    out = runs / 'pices_p7_accuracy'
    (out / 'cases').mkdir(parents=True, exist_ok=True)
    rows = []
    for n in (5, 9, 17):
        for method in ('pg', 'pices'):
            rows.append((n, method, 64, 32))
    rows += [(9, 'pices', p, 32) for p in (32, 128)]
    rows += [(17, m, 64, 64) for m in ('pg', 'pices')]
    manifest = []
    for n, method, particles, steps in rows:
        name = f'assim_{method}_n{n}_p{particles if method == "pices" else 0}_s{steps}'
        text = (runs / 'pices_p7/cases' / f'assim_{method}_dt_eighth.cfg').read_text()
        text = edit(text, dict(nodex=n, nodey=n, nodez=n, mgunitx=n-1, mgunity=n-1, mgunitz=n-1,
                              tracers_per_element=particles, maxstep=steps, maxtotstep=steps+1,
                              fixed_timestep=format(4e-7/steps, '.17g'), storage_spacing=steps,
                              CBF_frequency=steps, vlowstep=2000, perturblayer=(n+1)//2,
                              datafile='PICES_P7_ACCURACY'))
        (out / 'cases' / (name + '.cfg')).write_text(text)
        manifest.append(dict(name=name, scenario='assim', method=method, variant=name,
                             nodes=n, particles=particles if method == 'pices' else 0,
                             steps=steps, dt=4e-7/steps, prescribed=False, length_scale=1))
    for n in (5, 9, 17):
        (out / f'refstate_assim_{n}.txt').write_text(''.join(
            f'{1+.1*(n-1-i)/(n-1):.17g} 1 {.8-.4*i/(n-1):.17g} 1 {1+.2*(n-1-i)/(n-1):.17g}\n'
            for i in range(n)))
        forcing = out / f'forcing_{n}'
        forcing.mkdir(exist_ok=True)
        for age in range(4):
            (forcing / f'trench.{age}.xyz').write_text('')
            for cap in range(12):
                (forcing / f'age.{age}.{cap}').write_text((str(30+10*age)+'\n')*(n*n))
    (out / 'matrix.json').write_text(json.dumps(dict(schema=1, stage='P7_accuracy', cases=manifest,
        production_switch=False, final_time=4e-7,
        interpretation='PG is a comparator, not truth. Review within-method spatial, particle and time sensitivity together.'), indent=2)+'\n')
    (out / 'matrix.tsv').write_text(''.join(f'{r["name"]}\t{r["nodes"]}\t{r["steps"]}\tassim\n' for r in manifest))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=Path)
    print(prepare(parser.parse_args().runs))
