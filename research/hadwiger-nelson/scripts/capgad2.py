"""Extract the subgraph that actually caps de Grey's ring.

Sa is 397 points, but the statement that matters -- its D = 4 ring takes at
most two colours in every 4-colouring -- surely does not need all of them.
Whatever does need is the real gadget of the construction and the thing worth
rebuilding one level up; the rest is scaffolding the twelve-fold closure
brought along.

Peeling needs the test to be cheap, and rebuilt from scratch it is not: the
cap is an UNSAT proof and takes a minute on the full graph.  So the graph is
built ONCE with a selector per vertex -- a vertex whose selector is false is
required to carry no colour at all, which removes it and every constraint it
imposes -- and each trial is a single incremental call under assumptions.  The
solver keeps everything it learned between trials, which is the whole point.

Sweeping a fixed snapshot by identity and repeating until a sweep changes
nothing makes the result irreducible for the order, not merely order-dependent.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 2
t0 = time.time()
P = build_Sa(K)
n = len(P)
ZERO = Point(K.zero(), K.zero())
b = IntBasis.covering(P)
r = b.rows(P)
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
zi = P.index(ZERO)
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
print(f"Sa: {n} pts, {len(E)} edges, ring {ring}  [{time.time()-t0:.0f}s]",
      flush=True)

nv = n * k
sel = [nv + 1 + v for v in range(n)]
ind = [nv + n + 1 + c for c in range(k)]
top = nv + n + k
cls = []
for v in range(n):
    cls.append([-sel[v]] + [1 + v * k + c for c in range(k)])
    for c in range(k):
        cls.append([sel[v], -(1 + v * k + c)])
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
for c in range(k):
    cls.append([-ind[c]] + [1 + v * k + c for v in ring])
enc = CardEnc.atleast(lits=list(ind), bound=CAP + 1, top_id=top,
                      encoding=EncType.seqcounter)
sv = Solver(name="cd15", bootstrap_with=cls + list(enc.clauses))
print(f"solver built, {len(cls)+len(enc.clauses)} clauses"
      f"  [{time.time()-t0:.0f}s]", flush=True)

pinned = {zi} | set(ring)
present = set(range(n))


def holds(pres):
    a = [sel[v] if v in pres else -sel[v] for v in range(n)]
    return not sv.solve(assumptions=a)


assert holds(present), "the cap does not hold on Sa itself"
print(f"cap confirmed on the full graph  [{time.time()-t0:.0f}s]", flush=True)

rounds = 0
while True:
    rounds += 1
    snapshot = sorted(present - pinned)
    dropped = 0
    for v in snapshot:
        if v not in present:
            continue
        if holds(present - {v}):
            present.discard(v)
            dropped += 1
    print(f"  sweep {rounds}: dropped {dropped}, now {len(present)} points"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not dropped:
        break

keep = sorted(present)
kept_edges = [(a, c) for a, c in E if a in present and c in present]
print(f"\nirreducible: {len(keep)} points, {len(kept_edges)} edges"
      f"  [{time.time()-t0:.0f}s]", flush=True)
json.dump({"k": k, "cap": CAP,
           "points": [[[str(x) for x in P[i].x.c],
                       [str(y) for y in P[i].y.c]] for i in keep],
           "ring_positions": [keep.index(i) for i in ring],
           "centre_position": keep.index(zi)},
          open("capgadget.json", "w"))
print("written capgadget.json", flush=True)
print("DONE", flush=True)
