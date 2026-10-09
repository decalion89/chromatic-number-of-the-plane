"""Referee check 9: clause-by-clause validation of the CNF written by check_small.py, and its meaning.

Every clause must be one of
  (A) x(v,0) v x(v,1) v x(v,2)                       for each vertex v (at least one colour)
  (B) -x(i,k) v -x(j,k)                              for each edge ij and colour k (proper)
  (T) -x(a,s) v -x(b,s+1) v t(a,b)                   for each arc (a,b) of a listed cycle, all s
  (C) -t(a0,a1) v ... v -t(a_{m-1},a0)               for each listed cycle
with x(v,k) = 3v + k + 1, and every such clause must be present.  Soundness: from a model pick for each vertex
any true colour; (B) makes it proper, (T) makes t(a,b) true whenever that arc is tight, and (C) leaves a false t,
hence a non-tight arc, on every listed cycle.  So UNSAT => every proper 3-colouring makes a listed cycle tight.
Also: dropping each cycle clause in turn, is the formula still unsatisfiable (is every listed cycle needed)?
"""
import gzip, json, sys, itertools
from pysat.solvers import Solver

BASE, CNFDIR = sys.argv[1], sys.argv[2]


def load(name):
    with gzip.open(f'{BASE}/{name}', 'rt') as f:
        return json.load(f)


for wname, cnf in [('witness_q7.json.gz', 'q7.cnf'), ('witness_q31.json.gz', 'q31.cnf')]:
    W = load(wname); n = len(W['points']); E = [tuple(e) for e in W['edges']]; cycles = W['cycles']
    lines = open(f'{CNFDIR}/{cnf}').read().split('\n')
    hdr = lines[0].split(); nv, nc = int(hdr[2]), int(hdr[3])
    clauses = []
    for ln in lines[1:]:
        if ln.strip():
            lits = list(map(int, ln.split()))
            assert lits[-1] == 0
            clauses.append(tuple(lits[:-1]))
    x = lambda v, k: 3 * v + k + 1
    # identify the tight variables from the (T) clauses
    arcs = sorted({(cy[k], cy[(k + 1) % len(cy)]) for cy in cycles for k in range(len(cy))})
    tvar = {}
    expected = set()
    for v in range(n):
        expected.add(tuple(sorted([x(v, k) for k in range(3)])))
    for i, j in E:
        for k in range(3):
            expected.add(tuple(sorted([-x(i, k), -x(j, k)])))
    # find tight vars: a (T) clause has two negative colour literals and one positive literal > 3n
    for cl in clauses:
        pos = [l for l in cl if l > 3 * n]
        if len(cl) == 3 and len(pos) == 1:
            negs = sorted(-l for l in cl if l < 0)
            (a, ka), (b, kb) = [((u - 1) // 3, (u - 1) % 3) for u in negs]
            # orientation: tail a with colour s, head b with colour s+1
            if (kb - ka) % 3 == 1:
                arc = (a, b)
            else:
                arc = (b, a)
            tvar.setdefault(arc, set()).add(pos[0])
    ok_t = all(len(s) == 1 for s in tvar.values()) and len({next(iter(s)) for s in tvar.values()}) == len(tvar)
    tv = {a: next(iter(s)) for a, s in tvar.items()}
    for (a, b) in arcs:
        for s in range(3):
            expected.add(tuple(sorted([-x(a, s), -x(b, (s + 1) % 3), tv.get((a, b), 0)])))
    for cy in cycles:
        expected.add(tuple(sorted(-tv.get((cy[k], cy[(k + 1) % len(cy)]), 0) for k in range(len(cy)))))
    got = {tuple(sorted(c)) for c in clauses}
    print(f'== {cnf}: header {nv} vars {nc} clauses; parsed {len(clauses)} clauses ({len(got)} distinct)')
    print('   tight variables one per arc, distinct:', ok_t, '; arcs with a tight variable == arcs of listed cycles:',
          sorted(tv) == arcs, f'({len(arcs)} arcs)')
    print('   variables used = 3n + #arcs:', nv == 3 * n + len(arcs),
          '; clause set == expected set exactly:', got == expected, '; no duplicates:', len(got) == len(clauses))
    # semantic check: SAT without the cycle clauses, UNSAT with all
    cyc_cls = [tuple(sorted(-tv[(cy[k], cy[(k + 1) % len(cy)])] for k in range(len(cy)))) for cy in cycles]
    base = [c for c in clauses if tuple(sorted(c)) not in set(cyc_cls)]
    with Solver(name='cadical153', bootstrap_with=base) as s:
        print('   without cycle clauses: SAT =', s.solve())
    with Solver(name='glucose4', bootstrap_with=clauses) as s:
        print('   full formula (glucose4): SAT =', s.solve())
    need = []
    for drop in range(len(cycles)):
        cl2 = base + [c for k, c in enumerate(cyc_cls) if k != drop]
        with Solver(name='cadical153', bootstrap_with=cl2) as s:
            need.append(s.solve())
    print('   dropping listed cycle k makes it SAT (k = 0..):', need)
