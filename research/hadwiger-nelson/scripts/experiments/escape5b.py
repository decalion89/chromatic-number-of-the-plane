"""Iterate the five-colour cover: is one pair really enough?

The sampled hitting set said a single pair, (281, 655) at squared distance 3,
meets all 32 escapes of a 12-point neighbourhood in the 803-point graph.  A
sample only proves NECESSITY.  Iterating closes it: assume the cover's pairs
differ -- as clauses -- and ask for another escape.  When none survives, the
cover is SUFFICIENT and the requirement is exact.

Two ways it can end badly, both worth knowing:

  * the graph plus the cover stops being 5-colourable.  Then no gadget forcing
    those pairs apart can exist alongside this graph -- which, read the other
    way, says the graph plus such a gadget already refuses five.
  * the cover keeps growing.  Then the requirement is diffuse and no single
    gadget will carry it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
K = 5
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
NB = json.load(open(f"{ROOT}/data/escape5_five_247_c.json"))["neighbourhood"]
print(f"  n={n}, |N(p)|={len(NB)}   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            base.append([-X(v, a), -X(v, b)])
for x, y in g.edges():
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])
NBset = set(NB)
adj = {(a, b) for a in NB for b in g.adj[a] if b in NBset}
cover, rounds, total = [], 0, 0
while rounds < 30:
    rounds += 1
    cnf = list(base)
    for (a, b) in cover:
        for c in range(K):
            cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if not s.solve():
        print("  round %d: graph + cover is NOT 5-colourable -- a gadget "
              "forcing these %d pairs apart would refuse five by itself"
              % (rounds, len(cover)), flush=True)
        for (a, b) in cover:
            dd = (g.vertices[a] - g.vertices[b]).norm2()
            print("     %d,%d  d^2 = %s  (%.4f)" % (a, b, dd, float(dd)),
                  flush=True)
        s.delete(); break
    got = []
    status = None
    for _ in range(80):
        s.conf_budget(4_000_000)
        status = s.solve_limited(assumptions=[-X(u, c) for u in NB
                                              for c in (0, 1, 2)])
        if status is not True:
            break
        pos = set(l for l in s.get_model() if l > 0)
        got.append({u: next(c for c in range(K) if X(u, c) in pos) for u in NB})
        s.add_clause([-X(u, got[-1][u]) for u in NB])
    s.delete()
    total += len(got)
    if not got and status is None:
        # solve_limited returns None when the conflict budget runs out:
        # that decides nothing, so it must not be reported as "no escape".
        print("  round %d: UNDECIDED -- the solver ran out of its conflict "
              "budget before finding an escape or proving that none exists."
              % rounds, flush=True)
        break
    if not got:
        print("", flush=True)
        print("  round %d: NO ESCAPE SURVIVES.  Forcing these %d pair(s) apart "
              "raises mu_5 above 2." % (rounds, len(cover)), flush=True)
        for (a, b) in cover:
            dd = (g.vertices[a] - g.vertices[b]).norm2()
            print("     %d,%d  d^2 = %s  (%.4f)" % (a, b, dd, float(dd)),
                  flush=True)
        print("  %d escapes seen in total   [%.0fs]" % (total, time.time()-t0),
              flush=True)
        json.dump({"graph": "five_247_c.json", "neighbourhood": NB,
                   "rounds": rounds, "escapes_seen": total,
                   "cover": [[a, b, str((g.vertices[a]-g.vertices[b]).norm2())]
                             for a, b in cover]},
                  open(f"{ROOT}/data/mu5_requirement.json", "w"))
        print("  written data/mu5_requirement.json", flush=True)
        break
    sets = []
    for col in got:
        mono = {(a, b) for i, a in enumerate(NB) for b in NB[i+1:]
                if col[a] == col[b] and (a, b) not in adj and (b, a) not in adj}
        sets.append(mono)
    if any(not m for m in sets):
        print("  round %d: an escape has no monochromatic non-adjacent pair -- "
              "pairs alone cannot raise mu_5 here" % rounds, flush=True)
        break
    unc = list(range(len(sets)))
    added = 0
    while unc:
        cnt = Counter()
        for i in unc:
            for pr in sets[i]:
                cnt[pr] += 1
        pr, _ = cnt.most_common(1)[0]
        cover.append(pr); added += 1
        unc = [i for i in unc if pr not in sets[i]]
    print("  round %d: %d escapes, +%d pairs, cover now %d   [%.0fs]"
          % (rounds, len(got), added, len(cover), time.time()-t0), flush=True)
