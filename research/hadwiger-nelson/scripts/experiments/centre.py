"""The centre of a unit triangle, where 2 + 3 = 5.

A complete reduction, and the argument is four lines.  Let h be a vertex, T a
unit triangle whose centre is h -- so its three vertices sit at distance
1/sqrt3 from h, the circumradius of a unit triangle.  Suppose

    (P)   every 5-colouring gives  c(h) = c(t)  for some t in T.

Let rho be the rotation by 120 degrees about h and take G* = G u rho G u rho^2 G.
Rho fixes h, so each copy satisfies P about the SAME h: for each i there is
s_i with c(rho^i t_{s_i}) = c(h).  Three copies, three indices, so two copies
i != j share an index s -- and rho^i t_s, rho^j t_s are two points of the same
unit triangle, hence adjacent, hence cannot both be c(h).  Contradiction.

    (P) at any vertex of any unit-distance graph  =>  chi(R^2) >= 6,
    and the witness is three rotated copies.

Nothing about T's position is used beyond rho permuting it, so this needs no
tuned angle, no field extension and no search over rotations.

What makes (P) reachable is that its two halves are each AT their local ceiling
and add to exactly five.  c(h) avoids c(N(h)); (P) says it also avoids c(T);
so (P) fails only when some colour escapes both.  T is a triangle, so it always
shows 3 -- free.  N(h) lies on the unit circle, where chi is 2 and no graph can
do better.  3 + 2 = 5, exactly, with nothing to spare:

    mu_5( N(h) u T ) = 5   is equivalent to (P), and is what to measure.

Every mu in this project was measured on N(h) alone, where the ceiling is 2 and
the gap to five is three rungs.  On N(h) u T the triangle pays three of them.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAMES = sys.argv[1:] or ["five_247_c.json", "five_247.json"]
for NAME in NAMES:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    E = list(g.edges())
    adj = defaultdict(set)
    for x, y in E:
        adj[x].add(y); adj[y].add(x)
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    one = F.rational(Fr(1)); third = F.rational(Fr(1, 3))
    ring = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (g.vertices[i]-g.vertices[j]).norm2() == third:
                ring[i].append(j); ring[j].append(i)
    # the triangles centred at h: three ring points pairwise a unit apart
    tris = defaultdict(list)
    for h, vs in ring.items():
        for a in range(len(vs)):
            for b in range(a + 1, len(vs)):
                if (g.vertices[vs[a]] - g.vertices[vs[b]]).norm2() != one: continue
                for c in range(b + 1, len(vs)):
                    if (g.vertices[vs[a]] - g.vertices[vs[c]]).norm2() == one and \
                       (g.vertices[vs[b]] - g.vertices[vs[c]]).norm2() == one:
                        tris[h].append((vs[a], vs[b], vs[c]))
    print(f"\n{NAME} n={n}: {len(tris)} vertices centre a unit triangle; "
          f"triangle counts {Counter(len(v) for v in tris.values())}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not tris: continue
    X = lambda v, c: 1 + v * K + c
    base = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in E:
        for c in range(K):
            base.append([-X(x, c), -X(y, c)])
    ms = MuSolver(g, 5, budget=None)
    print(f"  base colouring {ms.colourable}   [{time.time()-t0:.0f}s]", flush=True)
    seen = Counter(); best = None
    order = sorted(tris, key=lambda h: -(len(adj[h]) + 3 * len(tris[h])))
    for h in order[:150]:
        for T in tris[h][:3]:
            XT = sorted(set(adj[h]) | set(T))
            v = ms.mu(XT)
            seen[v] += 1
            if best is None or (v or 0) > best[0]:
                best = (v, h, len(adj[h]), len(XT))
                print(f"    mu(N(h) u T) = {v}  at h={h}, |N(h)|={len(adj[h])}, "
                      f"|X|={len(XT)}   [{time.time()-t0:.0f}s]", flush=True)
            if v == 5:
                print(f"    *** (P) HOLDS at h={h} with T={T} ***", flush=True)
                json.dump({"graph": NAME, "h": h, "T": list(T)},
                          open(f"{ROOT}/data/centre_hit.json", "w"))
                sys.exit(0)
    print(f"  mu(N(h) u T) over {sum(seen.values())} pairs: "
          f"{dict(sorted(seen.items(), key=str))}   [{time.time()-t0:.0f}s]",
          flush=True)
    ms.close()
