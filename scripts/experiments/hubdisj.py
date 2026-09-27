"""A HUB disjunction, derived from a target point instead of hunted blind.

mu_4(G, p) = 3 says two things: every 4-colouring already puts three colours on
N(p), and at least one puts only three.  Those are the ESCAPES.  To block p --
to make G + p refuse four -- every escape has to be destroyed.

An escape is a proper colouring whose restriction to N(p) uses three colours,
so it partitions N(p) into three classes, and any two members of the same class
form a monochromatic pair.  Forcing either of those two apart destroys that
escape.  Two points of N(p) that are already adjacent (60 degrees apart on the
circle) can never be monochromatic and are never candidates.

So the question "which points must be forced?" is a HITTING SET problem, stated
exactly:

    collect the escapes;
    each contributes the set of its monochromatic non-adjacent pairs in N(p);
    a family of pairs blocks p only if it meets every one of those sets.

A hitting set over a SAMPLE of escapes is a necessary condition, not a
sufficient one -- destroying the sampled escapes is required of any solution,
and more escapes may remain.  That is the useful direction: it gives a lower
bound on what a gadget would have to do, and it names the pairs, and with them
their distances, which is what decides whether such a gadget can exist at all.

The frequency of each pair across escapes is the other half of the answer: a
pair that is monochromatic in almost every escape is load-bearing, and a pair
that is monochromatic in one is a detail.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
F = Field((3, 11, 247))
one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
Sa = build_Sa(F); g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
seen, U = set(Sa), list(Sa)
for w in orbit(Sa[265]):
    rot = g60.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
print(f"carrier n={n} m={sum(len(a) for a in g.adj)//2} deg="
      f"{2.0*sum(len(a) for a in g.adj)/2/n:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
S = set(g.vertices)
cells = defaultdict(list)
for i, p in enumerate(g.vertices):
    cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
best = None
for c in deg_order[:40]:
    for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
        for p in g.vertices:
            z = rot(p)
            if z in S:
                continue
            zx, zy = z.fx, z.fy
            cx, cy = int(zx // 1), int(zy // 1)
            nb = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for i in cells.get((cx + dx, cy + dy), ()):
                        ex, ey = gx[i] - zx, gy[i] - zy
                        if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                           (g.vertices[i] - z).norm2() == one:
                            nb.append(i)
            if best is None or len(nb) > len(best[1]):
                best = (z, sorted(nb))
    if best and len(best[1]) >= 15:
        break
P, NB = best
print(f"  target |N(p)| = {len(NB)}   [{time.time()-t0:.0f}s]", flush=True)

K = 4
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            base.append([-X(v, a), -X(v, b)])
for x, y in g.edges():
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

# The four-pair disjunction the cover produced has no shared endpoint, and an
# earlier pass already knew why that is fatal: each pair carries its own centre
# of rotation and one rotation cannot serve them all.  The stack only consumes
# a HUB disjunction --
#
#     in every k-colouring,  c(h) = c(q)  for some q in W
#
# -- because then the rotation about h repeats the same statement in every
# copy, over a shifted W.  So ask for that shape directly.  For each h in
# N(p), every escape contributes the set of members sharing h's colour, and a
# hitting set over those sets is exactly a W that works.  An escape where h's
# colour is unique in N(p) kills that hub outright: no W can catch it.
s = Solver(name="cd19", bootstrap_with=base)
assert s.solve()
escapes = []
for _ in range(200):
    s.conf_budget(4_000_000)
    if s.solve_limited(assumptions=[-X(u, 0) for u in NB]) is not True:
        break
    pos = set(l for l in s.get_model() if l > 0)
    escapes.append({u: next(c for c in range(K) if X(u, c) in pos) for u in NB})
    s.add_clause([-X(u, escapes[-1][u]) for u in NB])
s.delete()
print("  %d escapes sampled   [%.0fs]" % (len(escapes), time.time() - t0),
      flush=True)

best = None
for h in NB:
    sets = []
    dead = False
    for col in escapes:
        same = {q for q in NB if q != h and col[q] == col[h]}
        if not same:
            dead = True
            break
        sets.append(same)
    if dead:
        continue
    unc = list(range(len(sets)))
    W = []
    while unc:
        cnt = Counter()
        for i in unc:
            for q in sets[i]:
                cnt[q] += 1
        q, _ = cnt.most_common(1)[0]
        W.append(q)
        unc = [i for i in unc if q not in sets[i]]
    if best is None or len(W) < len(best[1]):
        best = (h, W)
    print("    hub %d: W of size %d" % (h, len(W)), flush=True)
if best is None:
    print("  every hub is killed by some escape: no hub disjunction here",
          flush=True)
else:
    h, W = best
    print("", flush=True)
    print("  smallest hub disjunction: c(%d) meets one of %d points"
          % (h, len(W)), flush=True)
    for q in W:
        d2 = (g.vertices[h] - g.vertices[q]).norm2()
        print("     %d at d^2 = %s   (%.4f)" % (q, d2, float(d2)), flush=True)
    rr = {str((g.vertices[h] - g.vertices[q]).norm2()) for q in W}
    print("  distinct radii from the hub: %d %s" % (len(rr), sorted(rr)),
          flush=True)
    json.dump({"n": n, "hub": h, "W": list(W),
               "radii": [str((g.vertices[h] - g.vertices[q]).norm2())
                         for q in W]},
              open(ROOT + "/data/hub_disjunction.json", "w"))
    print("  written data/hub_disjunction.json", flush=True)
