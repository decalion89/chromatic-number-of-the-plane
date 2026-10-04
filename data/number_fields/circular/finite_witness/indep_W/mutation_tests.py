#!/usr/bin/env python3
"""Negative controls for the referee checkers: each mutation of witness_q455 must be rejected.
Usage: python3 mutation_tests.py FW_DIR OUT_DIR"""
import gzip, json, os, subprocess, sys, copy

FW, OUT = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
W0 = json.load(gzip.open(os.path.join(FW, 'witness_q455.json.gz'), 'rt'))
C0 = json.load(gzip.open(os.path.join(FW, 'witness_q455_critical.json.gz'), 'rt'))


def dump(obj, name):
    p = os.path.join(OUT, name)
    with gzip.open(p, 'wt') as f:
        json.dump(obj, f)
    return p


def run(args):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr).strip().splitlines()[-1][:160]


muts = {}
W = copy.deepcopy(W0); W['points'][5][0] += 1; muts['point coordinate +1'] = W
W = copy.deepcopy(W0); W['edges'].pop(17); muts['edge removed (not induced)'] = W
W = copy.deepcopy(W0)
E = {tuple(sorted(e)) for e in W['edges']}
pair = next((i, j) for i in range(len(W['points'])) for j in range(i + 1, len(W['points'])) if (i, j) not in E)
W['edges'].append(list(pair)); muts['non-unit pair added as edge'] = W
W = copy.deepcopy(W0); i, j = W['edges'][0]; W['colouring'][j] = (W['colouring'][i] + 1) % 7; muts['colour changed'] = W
W = copy.deepcopy(W0); C = W['cycles'][0]; C[1], C[2] = C[2], C[1]; muts['cycle vertices swapped'] = W
W = copy.deepcopy(W0); W['points'][3] = list(W['points'][4]); muts['duplicate point'] = W
W = copy.deepcopy(W0); W['d'] = 4 * 455; muts['d not squarefree'] = W
ok = True
for name, W in muts.items():
    p = dump(W, 'mut.json.gz')
    rc, last = run([os.path.join(HERE, 'referee_check.py'), p, os.path.join(OUT, 'cnf')])
    print(f"referee_check  | {name:32s} | rejected={rc != 0} | {last}")
    ok &= rc != 0
# criticality certificates
wp = dump(W0, 'w.json.gz')
# (a) replace c_v by the stored colouring of H (which has a tight cycle) for a vertex v off that cycle
col = W0['colouring']
tight = [C for C in W0['cycles'] if all((col[C[(k+1) % len(C)]] - col[C[k]]) % 7 == 2 for k in range(len(C)))][0]
v = next(u for u in range(len(col)) if u not in tight)
C = copy.deepcopy(C0); C['critical_colourings'][v] = [(-1 if u == v else col[u]) for u in range(len(col))]
rc, last = run([os.path.join(HERE, 'referee_critical.py'), wp, dump(C, 'c.json.gz')])
print(f"referee_critical | stored colouring as c_{v} (tight cycle) | rejected={rc != 0} | {last}"); ok &= rc != 0
C = copy.deepcopy(C0); c = C['critical_colourings'][0]; e = next(e for e in W0['edges'] if 0 not in e)
c[e[1]] = c[e[0]]; rc, last = run([os.path.join(HERE, 'referee_critical.py'), wp, dump(C, 'c.json.gz')])
print(f"referee_critical | improper c_0 | rejected={rc != 0} | {last}"); ok &= rc != 0
C = copy.deepcopy(C0); C['critical_colourings'][3][3] = 0
rc, last = run([os.path.join(HERE, 'referee_critical.py'), wp, dump(C, 'c.json.gz')])
print(f"referee_critical | c_3(3) != -1 | rejected={rc != 0} | {last}"); ok &= rc != 0
C = copy.deepcopy(C0); C['critical_colourings'].pop()
rc, last = run([os.path.join(HERE, 'referee_critical.py'), wp, dump(C, 'c.json.gz')])
print(f"referee_critical | one certificate missing | rejected={rc != 0} | {last}"); ok &= rc != 0
print("ALL MUTATIONS REJECTED" if ok else "SOME MUTATION ACCEPTED")
