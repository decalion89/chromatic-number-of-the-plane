"""Shrink G by unsatisfiable core, then by hand.

G is 1581 vertices and not 4-colourable, which is de Grey's theorem.  Most of
those vertices are scaffolding: the property is an UNSAT proof, and an UNSAT
proof under assumptions hands back the subset of assumptions it actually used.
Giving every vertex a selector -- false means the vertex carries no colour at
all, which removes it and every constraint it imposes -- turns "which vertices
does the proof need" into one solver call instead of 1581.

Core-guided passes first, until the core stops shrinking, then one vertex at a
time until a whole sweep changes nothing.  Every intermediate is verified the
same way it was produced: the graph is rebuilt from the surviving points with
exact arithmetic and re-solved, so a bug in the selector encoding cannot
survive into the answer.

The smallest 5-chromatic unit-distance graph known is 509 vertices, so the
number this lands on is worth having either way.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
P = build_G(K, as_graph=False)
n = len(P)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
print(f"G: {n} pts, {len(E)} edges, headroom {b.overflow_headroom(r):.2e}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

nv = n * k
sel = [nv + 1 + v for v in range(n)]
cls = []
for v in range(n):
    cls.append([-sel[v]] + [1 + v * k + c for c in range(k)])
    for c in range(k):
        cls.append([sel[v], -(1 + v * k + c)])
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
sv = Solver(name="cd15", bootstrap_with=cls)
print(f"solver built, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)


def blocks(pres, want_core=False):
    a = [sel[v] if v in pres else -sel[v] for v in range(n)]
    if sv.solve(assumptions=a):
        return (False, None) if want_core else False
    if not want_core:
        return True
    core = set(sv.get_core() or [])
    return True, {v for v in range(n) if sel[v] in core}


def verify(pres):
    """Rebuild from the surviving points and re-solve, trusting nothing."""
    Q = [P[i] for i in sorted(pres)]
    bb = IntBasis.covering(Q)
    rr = bb.rows(Q)
    assert bb.overflow_headroom(rr) < 1.0
    EE = sorted(set((min(a, c), max(a, c))
                    for a, c in fast_edges_complete(bb, rr)))
    m = len(Q)
    cc = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, c in EE:
        for col in range(k):
            cc.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cc)
    ok = s.solve()
    s.delete()
    return (not ok), m, len(EE)


ok, core = blocks(set(range(n)), want_core=True)
assert ok, "G is 4-colourable?!"
print(f"first core: {len(core)} of {n} vertices  [{time.time()-t0:.0f}s]",
      flush=True)
present = set(core)
while True:
    ok, c2 = blocks(present, want_core=True)
    assert ok, "the core stopped blocking -- encoding bug"
    print(f"  core pass: {len(present)} -> {len(c2)}  [{time.time()-t0:.0f}s]",
          flush=True)
    if len(c2) >= len(present):
        break
    present = set(c2)
good, m, ne = verify(present)
print(f"after cores: {m} points, {ne} edges, still blocks: {good}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
assert good

rounds = 0
while True:
    rounds += 1
    snapshot = sorted(present)
    dropped = 0
    for v in snapshot:
        if v not in present:
            continue
        if blocks(present - {v}):
            present.discard(v)
            dropped += 1
    good, m, ne = verify(present)
    print(f"  sweep {rounds}: dropped {dropped}, now {m} points, {ne} edges, "
          f"blocks={good}  [{time.time()-t0:.0f}s]", flush=True)
    assert good
    json.dump({"k": k, "points": [[[str(x) for x in P[i].x.c],
                                   [str(y) for y in P[i].y.c]]
                                  for i in sorted(present)]},
              open("ming.json", "w"))
    if not dropped:
        break
print(f"\nirreducible: {m} points, {ne} edges  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
