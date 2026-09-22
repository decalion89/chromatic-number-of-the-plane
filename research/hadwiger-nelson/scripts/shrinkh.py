"""Shrink the ingredient, not the result.

Z refuses four colours because H forces H[157] and H[327] to share a colour
and the rotation puts the two images of H[327] one apart.  So any subgraph
H' of H that still forces that pair gives a smaller Z' = H' u rho(H') that
still refuses four colours.  Shrinking H' is the easier problem by far: the
instance is "H' is 4-colourable with c(a)=0 and c(b)=1", which is heavily
constrained and proves UNSAT quickly, where "Z' is 4-colourable" is not.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

K1 = Field((3, 11, 247))
pts = build_Sa(K1)
r1 = rotation_joining(Fr(1), K1).about(pts[25])
seen, H = set(), []
for p in pts:
    for q in (p, r1(p)):
        if q not in seen: seen.add(q); H.append(q)
g = build_graph(H)
n, K = g.n, 4
A_, B_ = 157, 327
print(f"H: n={n} m={sum(len(a) for a in g.adj)//2}; forcing pair {A_},{B_}", flush=True)
X = lambda v, c: 1 + v * K + c
SEL = lambda v: 1 + n * K + v
cnf = [[-SEL(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
cnf.append([X(A_, 0)])
cnf.append([X(B_, 1)])

def forces(sub):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[SEL(v) for v in sorted(sub)])
    if ok:
        s.delete(); return None
    core = s.get_core(); s.delete()
    return set(l - 1 - n * K for l in core) if core else set(sub)

t0 = time.time()
sub = set(range(n))
assert forces(sub) is not None
for r in range(80):
    c = forces(sub | {A_, B_})
    if c is None or len(c) >= len(sub): break
    sub = c | {A_, B_}
    print(f"  core round {r+1}: {len(sub)}   [{time.time()-t0:.0f}s]", flush=True)
print(f"  after cores: {len(sub)}   [{time.time()-t0:.0f}s]", flush=True)
order = sorted(sub, key=lambda v: len(g.adj[v]))
dropped = 0
for v in order:
    if v in (A_, B_) or v not in sub: continue
    trial = sub - {v}
    c = forces(trial)
    if c is not None:
        sub = (c | {A_, B_}) if len(c | {A_, B_}) < len(trial) else trial
        dropped += 1
print(f"\nminimal forcing subgraph H': {len(sub)} vertices "
      f"(dropped {dropped})   [{time.time()-t0:.0f}s]", flush=True)

Hp = [H[v] for v in sorted(sub)]
idx = {v: i for i, v in enumerate(sorted(sub))}
rot = rotation_joining(Fr(64, 9), K1).about(H[A_])
seen, Zp = set(), []
for p in Hp:
    for q in (p, rot(p)):
        if q not in seen: seen.add(q); Zp.append(q)
gz = build_graph(Zp)
mz = sum(len(a) for a in gz.adj) // 2
XX = lambda v, c: 1 + v * K + c
cc = [[XX(v, c) for c in range(K)] for v in range(gz.n)]
for u, v in gz.edges():
    for c in range(K):
        cc.append([-XX(u, c), -XX(v, c)])
s = Solver(name="m22", bootstrap_with=cc); four = s.solve(); s.delete()
print(f"Z' = H' u rho(H'): n={gz.n} m={mz}  4-colourable={four}", flush=True)
if not four:
    print(f"  *** a {gz.n}-vertex 5-chromatic unit-distance graph ***", flush=True)
    def dump(pt):
        return [[[c.numerator, c.denominator] for c in pt.x.c],
                [[c.numerator, c.denominator] for c in pt.y.c]]
    json.dump({"field_generators": list(K1.gens), "n": gz.n,
               "from": "H' u rho(H') with H' the minimal subgraph of H forcing the pair",
               "points": [dump(p) for p in gz.vertices]},
              open("/home/user/darwin-50/research/hadwiger-nelson/data/five_247_small.json", "w"))
    print("  written data/five_247_small.json", flush=True)
