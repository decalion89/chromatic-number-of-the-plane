"""The three-fold disjunctive spindle, aimed at five.

The gadget: if every 5-colouring of H has  c(v) = c(q1) or c(v) = c(q2)  with
|v - q1| = |v - q2| = 1/sqrt3, then H u rho(H) u rho^2(H) for rho the 120
degree rotation about v has NO 5-colouring, because three copies must each
pick one of two options, two of them pick the same, and the two corresponding
points are adjacent vertices of a unit equilateral triangle about v.

This is strictly weaker than the hypothesis de Grey's spindle needs -- a
forced pair is the special case q1 = q2 -- and it costs one extra copy.  It
needs no radical beyond sqrt3, since cos 120 = -1/2 and sin 120 = sqrt3/2.

Sa gave zero at four colours, exhaustively, over all 14 919 option-pairs.  The
question that matters is five, on the graphs that already refuse four.

Testing it is one warm SAT call per option-pair.  "c(v) differs from both" is
not a clause, but the CNF carries no symmetry breaking, so any colouring that
separates v from q1 and from q2 can be permuted to put c(v) = 0 -- and then the
condition is the two unit assumptions -X(q1,0), -X(q2,0).
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()

def load(name, gens):
    F = Field(gens)
    d = json.load(open(f"{ROOT}/data/{name}"))
    return F, [Point(F.element([Fr(a, b) for a, b in x]),
                     F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]

def hexsum(F, Z):
    half = F.rational(Fr(1, 2)); rt = F.sqrt(3) * half
    one, zero = F.rational(Fr(1)), F.rational(Fr(0))
    H = [Point(zero, zero), Point(one, zero), Point(half, rt), Point(-half, rt),
         Point(-one, zero), Point(-half, -rt), Point(half, -rt)]
    seen, out = set(), []
    for t in H:
        for p in Z:
            q = Point(p.x + t.x, p.y + t.y)
            if q not in seen: seen.add(q); out.append(q)
    return out

CASES = []
F1, Z = load("five_247_c.json", (3, 11, 247))
CASES.append(("five_247_c", F1, Z))
CASES.append(("five_247_c (+) H", F1, hexsum(F1, Z)))
F2, Zb = load("five_247.json", (3, 11, 247))
CASES.append(("five_247", F2, Zb))
F3, Zs = load("five_symmetric.json", (3, 11, 247))
CASES.append(("five_symmetric", F3, Zs))

for name, F, U in CASES:
    THIRD = F.rational(Fr(1, 3))
    g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
    at13 = defaultdict(list)
    for i in range(n):
        pi = g.vertices[i]
        for j in range(i + 1, n):
            if (pi - g.vertices[j]).norm2() == THIRD:
                at13[i].append(j); at13[j].append(i)
    hubs = {v: qs for v, qs in at13.items() if len(qs) >= 2}
    tot = sum(len(q) * (len(q) - 1) // 2 for q in hubs.values())
    print(f"  {name:20s} n={n} m={m} deg={2.0*m/n:.2f}: {len(hubs)} pivots "
          f"with two or more at d^2=1/3, {tot} option-pairs"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not hubs:
        continue
    K = 5
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for c in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, c)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if not s.solve():
        print(f"    *** not even {K}-colourable ***", flush=True)
        s.delete(); continue
    found, tried = [], 0
    for v, qs in hubs.items():
        for a in range(len(qs)):
            for b in range(a + 1, len(qs)):
                tried += 1
                if not s.solve(assumptions=[X(v, 0), -X(qs[a], 0), -X(qs[b], 0)]):
                    found.append((v, qs[a], qs[b]))
                    if len(found) <= 5:
                        print(f"    *** DISJUNCTIVE FORCING AT FIVE: "
                              f"v={v} q1={qs[a]} q2={qs[b]} ***", flush=True)
    s.delete()
    print(f"    tried {tried} option-pairs at five, {len(found)} disjunctions"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if found:
        v, q1, q2 = found[0]
        rho = Rotation(F.rational(Fr(-1, 2)), F.sqrt(3) * F.rational(Fr(1, 2)))
        rot = rho.about(g.vertices[v])
        seen, V, cur = set(g.vertices), list(g.vertices), list(g.vertices)
        for _ in range(2):
            cur = [rot(p) for p in cur]
            for p in cur:
                if p not in seen: seen.add(p); V.append(p)
        gg = build_graph(V); nn = gg.n; mm = sum(len(a) for a in gg.adj) // 2
        print(f"    three-fold union: n={nn} m={mm} deg={2.0*mm/nn:.2f}",
              flush=True)
        Y = lambda u, c: 1 + u * K + c
        cl = [[Y(u, c) for c in range(K)] for u in range(nn)]
        for u in range(nn):
            for a in range(K):
                for c in range(a + 1, K):
                    cl.append([-Y(u, a), -Y(u, c)])
        for x, y in gg.edges():
            for c in range(K):
                cl.append([-Y(x, c), -Y(y, c)])
        t1 = time.time()
        sv = Solver(name="cd19", bootstrap_with=cl); r = sv.solve(); sv.delete()
        print(f"    {K}-colourable: {r}   [{time.time()-t1:.0f}s]", flush=True)
        if not r:
            json.dump({"field_generators": list(F.gens), "n": nn, "m": mm,
                       "mechanism": "three-fold disjunctive spindle at 1/sqrt3",
                       "pivot_q1_q2": [v, q1, q2],
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in gg.vertices]},
                      open(f"{ROOT}/data/six_candidate.json", "w"))
            print("    *** written data/six_candidate.json ***", flush=True)
