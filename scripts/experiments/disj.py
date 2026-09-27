"""Disjunctive forcing, and the three-fold spindle that consumes it.

De Grey's spindle needs a pair FORCED equal: c(v) = c(q) in every k-colouring.
That is a very strong hypothesis, and the filter says our graphs do not have it
at five -- 1 320 vertices, 1.3 million pairs, every single one separated.

There is a weaker hypothesis that still closes the argument.  Suppose only that
in every k-colouring of H,

        c(v) = c(q1)   OR   c(v) = c(q2)

with |v - q1| = |v - q2| = r, and let rho be the rotation by 2*pi/3 about v.
Form  H u rho(H) u rho^2(H).  All three copies share v, which rho fixes, so in
each copy j the disjunction gives c(v) = c(rho^j(q_i)) for some i(j) in {1,2}.
Three copies, two options: two copies j != k take the SAME option i.  Then

        c(rho^j(q_i)) = c(v) = c(rho^k(q_i))

and rho^j(q_i), rho^k(q_i) are two vertices of an equilateral triangle
inscribed in the circle of radius r about v.  Its side is r*sqrt3, so if
r = 1/sqrt3 that side is exactly 1 and the two are ADJACENT.  Contradiction.
H u rho(H) u rho^2(H) refuses k colours.

Why three and why 1/sqrt3: with s copies and options at radius r_i the copies
that clash form a circulant on Z_s with connection set +-j_i, where the angle
2 arcsin(1/(2 r_i)) is 2 pi j_i / s.  An option's copies must be independent
there, so at most floor(s/2) copies per option, and 2*floor(s/2) < s exactly
when s is ODD.  s = 3 is the smallest, and it needs only cos 120 = -1/2 and
sin 120 = sqrt3/2 -- no radical the field does not already have.  (s = 5 also
works, at radii (5 +- sqrt5)/10, but sin 72 = sqrt(10 + 2 sqrt5)/4 is a nested
radical and the multiquadratic field cannot hold it.)

So the target drops from "a forced pair" to "a forced pair of options", and the
geometric price is one squared distance: 1/3.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa, build_G, build_Y
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
F = Field((3, 5, 7, 11))
THIRD = F.rational(Fr(1, 3))
rot60 = _rot60(F)
rho120 = Rotation(F.rational(Fr(-1, 2)), F.sqrt(3) * F.rational(Fr(1, 2)))

Sa = build_Sa(F)
def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out

def grow(base):
    seen, U = set(Sa), list(Sa)
    for w in orb(Sa[base]):
        rot = rotation_joining(Fr(1), F).about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    return U

CASES = [("Sa", list(Sa)), ("Y", list(build_Y(F))), ("G", list(build_G(F, as_graph=False)))]
for b in (1, 7, 25):
    CASES.append((f"Sa[{b}]", grow(b)))

for name, U in CASES:
    g = build_graph(U); n = g.n
    m = sum(len(a) for a in g.adj) // 2
    # every ordered pair at squared distance exactly 1/3
    at13 = defaultdict(list)
    for i in range(n):
        pi = g.vertices[i]
        for j in range(n):
            if i != j and (pi - g.vertices[j]).norm2() == THIRD:
                at13[i].append(j)
    hubs = {v: qs for v, qs in at13.items() if len(qs) >= 2}
    tot = sum(len(qs) * (len(qs) - 1) // 2 for qs in hubs.values())
    print(f"  {name:7s} n={n} m={m} deg={2.0*m/n:.2f}: "
          f"{len(at13)} vertices have a neighbour at d^2=1/3, "
          f"{len(hubs)} have two or more, {tot} candidate option-pairs"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not hubs:
        continue
    K = 4 if name in ("Sa", "Y") else 5
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
        print(f"    (not {K}-colourable at all -- skipping)", flush=True)
        s.delete(); continue
    # "c(v) differs from both" is not a clause, but colour symmetry turns it
    # into one line of assumptions: the CNF has no symmetry breaking, so any
    # colouring separating v from q1 and from q2 can be permuted to put
    # c(v) = 0, and the condition becomes  -X(q1,0) and -X(q2,0).  One warm
    # solver for all of them instead of fifteen thousand cold ones.
    found, tried = [], 0
    for v, qs in hubs.items():
        for a in range(len(qs)):
            for b2 in range(a + 1, len(qs)):
                q1, q2 = qs[a], qs[b2]
                tried += 1
                if not s.solve(assumptions=[X(v, 0), -X(q1, 0), -X(q2, 0)]):
                    found.append((v, q1, q2))
                    if len(found) <= 6:
                        print(f"    *** DISJUNCTIVE FORCING at {K}: "
                              f"v={v} q1={q1} q2={q2} ***", flush=True)
    s.delete()
    print(f"    tried {tried} option-pairs at {K} colours, "
          f"{len(found)} disjunctions   [{time.time()-t0:.0f}s]", flush=True)
    if found:
        v, q1, q2 = found[0]
        rot = rho120.about(g.vertices[v])
        seen, V = set(g.vertices), list(g.vertices)
        cur = list(g.vertices)
        for _ in range(2):
            cur = [rot(p) for p in cur]
            for p in cur:
                if p not in seen: seen.add(p); V.append(p)
        gg = build_graph(V); nn = gg.n
        mm = sum(len(a) for a in gg.adj) // 2
        print(f"    three-fold union: n={nn} m={mm} deg={2.0*mm/nn:.2f}",
              flush=True)
        Y2 = lambda v_, c: 1 + v_ * K + c
        cl = [[Y2(u, c) for c in range(K)] for u in range(nn)]
        for u in range(nn):
            for a in range(K):
                for c in range(a + 1, K):
                    cl.append([-Y2(u, a), -Y2(u, c)])
        for x, y in gg.edges():
            for c in range(K):
                cl.append([-Y2(x, c), -Y2(y, c)])
        t1 = time.time()
        sv = Solver(name="cd19", bootstrap_with=cl); r = sv.solve(); sv.delete()
        print(f"    {K}-colourable: {r}   [{time.time()-t1:.0f}s]", flush=True)
        if not r:
            tag = "six_candidate" if K == 5 else "five_threefold"
            json.dump({"field_generators": list(F.gens), "n": nn, "m": mm,
                       "mechanism": "three-fold disjunctive spindle",
                       "pivot_q1_q2": [v, q1, q2], "K_refused": K,
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in gg.vertices]},
                      open(f"{ROOT}/data/{tag}.json", "w"))
            print(f"    *** written data/{tag}.json ***", flush=True)
