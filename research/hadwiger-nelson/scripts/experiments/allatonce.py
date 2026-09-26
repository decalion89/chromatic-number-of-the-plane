"""Settle every candidate pair with one solver call instead of 393.

Asking "can u and v differ?" once per pair costs a hard SAT call each, and on
the 7141-vertex symmetric graph that is hours.  But the question has a much
better shape: if there is a single proper 5-colouring in which ALL the
candidate pairs differ, then every one of them differs in some colouring, so
none is forced -- and that is one call.

And it has a trivial encoding.  "u and v differ" is exactly the constraint an
edge imposes, so the question is just

    is  G + (all candidate pairs as edges)  still 5-colourable?

A yes is a complete certificate for the whole candidate set.  A no is not a
forced pair -- it only says they cannot all differ at once -- and then the set
is split and the halves asked separately, which is a bisection over a few
calls rather than a scan over hundreds.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

d = json.load(open(HN_DIR + "/data/five_symmetric.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
     for x, y in d["points"]]
g = build_graph(P); n = g.n; K = 5
cand = [tuple(x) for x in json.load(open(
    "/tmp/hn/cand7141.json"))]
print(f"n={n}, {len(cand)} candidate pairs", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        base.append([-X(u, c), -X(v, c)])
tri = g.find_clique(3)
for i, v in enumerate(tri):
    base.append([X(v, i)])
    for c in range(K):
        if c != i: base.append([-X(v, c)])

def separable(pairs):
    cnf = list(base)
    for a, b in pairs:
        for c in range(K):
            cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    ok = s.solve(); s.delete()
    return ok

# Bisecting down from the whole set pays the expensive direction first: a
# group of a thousand extra edges is a hard SAT instance, while UNSAT comes
# back instantly.  Starting from small chunks reverses that -- twenty-five
# extra edges are easy to satisfy -- and only a chunk that fails gets split.
CH = 25
t0 = time.time()
todo = [cand[i:i + CH] for i in range(0, len(cand), CH)][::-1]
proved, stuck = 0, []
while todo:
    grp = todo.pop()
    ok = separable(grp)
    if ok:
        proved += len(grp)
        if proved % 250 < CH:
            print(f"    proved separable: {proved}/{len(cand)}"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
        continue
    print(f"  a group of {len(grp)} cannot all differ; splitting"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if len(grp) == 1:
        a, b = grp[0]
        dd = (P[a] - P[b]).norm2()
        print(f"  *** v{a} v{b} CANNOT differ -- FORCED EQUAL, d^2={dd} ***",
              flush=True)
        stuck.append((a, b, str(dd)))
        continue
    h = len(grp) // 2
    todo.append(grp[:h]); todo.append(grp[h:])
print(f"\nproved separable: {proved} of {len(cand)};  forced: {len(stuck)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
json.dump(stuck, open("/tmp/hn/forced_allatonce.json", "w"))
