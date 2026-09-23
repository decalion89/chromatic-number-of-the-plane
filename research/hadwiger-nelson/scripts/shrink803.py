"""Degree-ordered batch removal: log-many solves per round, not one per vertex.

Greedy vertex-by-vertex minimisation needs one 4-colour UNSAT proof per vertex,
which is days at 803 vertices.  Batch removal gets the same monotone descent
with a handful of solves per round: sort by degree, try dropping the lowest
tenth in one go, and halve the batch until the rest still refuses four.  Each
accepted batch is a real cut; each rejected one costs a single solve and halves
the next attempt.

Low degree first because the object is a glue of a carrier with a rotated copy
plus one spindle: the forcing lives in the overlap, where the degrees are high,
and the periphery is decoration.  Nothing is moved, only deleted, so the graph
stays a unit-distance graph in the same field and the certificate survives.
"""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
PTS = [Point(F.element([Fr(a, b) for a, b in x]),
             F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
K = 4

def refuses_four(sub):
    g = build_graph([PTS[i] for i in sub]); n = g.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    # Pin a triangle.  Deciding colourability is invariant under permuting the
    # colours, so fixing one clique to 0,1,2 is free -- and it removes the 4!
    # symmetric copies of every refutation, which is what makes an UNSAT proof
    # at four expensive on a nearly-minimal graph.  (The union of 1320 points
    # refuted in 83 s; this one, with fewer constraints carrying the same
    # contradiction, is slower without the pin.)
    tri = g.find_clique(3) or []
    for i, x in enumerate(tri):
        cl.append([X(x, i)])
        for c in range(K):
            if c != i:
                cl.append([-X(x, c)])
    s = Solver(name="cd19", bootstrap_with=cl); r = s.solve(); s.delete()
    return (not r), g

cur = list(range(len(PTS)))
ok, g = refuses_four(cur)
print(f"start n={g.n} m={sum(len(a) for a in g.adj)//2} refuses4={ok}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
solves = 1
while True:
    g = build_graph([PTS[i] for i in cur])
    order = sorted(range(len(cur)), key=lambda i: len(g.adj[i]))
    batch = max(1, len(cur) // 10)
    cut = False
    while batch >= 1:
        drop = {cur[i] for i in order[:batch]}
        trial = [u for u in cur if u not in drop]
        ok, _ = refuses_four(trial); solves += 1
        if ok:
            print(f"  dropped {batch}: n={len(trial)}  ({solves} solves)"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
            cur = trial; cut = True
            gg = build_graph([PTS[i] for i in cur])
            json.dump({"field_generators": list(F.gens), "n": gg.n,
                       "m": sum(len(a) for a in gg.adj) // 2,
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in gg.vertices]},
                      open(f"{ROOT}/data/five_247_min.json", "w"))
            break
        batch //= 2
    if not cut:
        break
gg = build_graph([PTS[i] for i in cur])
mm = sum(len(a) for a in gg.adj) // 2
print(f"  FINAL n={gg.n} m={mm} deg={2.0*mm/gg.n:.2f}  after {solves} solves"
      f"   [{time.time()-t0:.0f}s]", flush=True)
