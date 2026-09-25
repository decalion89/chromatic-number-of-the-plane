"""The COMPLETE requirement: iterate the cover until no escape survives.

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
NBset = set(NB)
adj = {(a, b) for a in NB for b in g.adj[a] if b in NBset}

# A sampled hitting set is only NECESSARY: killing the escapes we happened to
# see says nothing about the ones we did not.  Iterating closes that gap.
# Assume the pairs found so far are forced apart -- as clauses, so the solver
# treats them as true -- and ask again for an escape.  Each new escape needs a
# new pair, the cover grows, and the loop ends when no escape survives.  At
# that point the cover is SUFFICIENT: forcing exactly those pairs apart blocks
# the point, and the run has proved it rather than sampled it.
cover = []
rounds = 0
total_escapes = 0
while rounds < 40:
    rounds += 1
    cnf = list(base)
    for (a, b) in cover:
        for c in range(K):
            cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    if not s.solve():
        print("  round %d: the carrier plus these constraints is not "
              "4-colourable at all -- the cover is too strong" % rounds,
              flush=True)
        s.delete(); break
    got = []
    status = None
    for _ in range(60):
        s.conf_budget(4_000_000)
        status = s.solve_limited(assumptions=[-X(u, 0) for u in NB])
        if status is not True:
            break
        pos = set(l for l in s.get_model() if l > 0)
        col = {u: next(c for c in range(K) if X(u, c) in pos) for u in NB}
        got.append(col)
        s.add_clause([-X(u, col[u]) for u in NB])
    s.delete()
    total_escapes += len(got)
    if not got and status is None:
        # solve_limited returns None when the conflict budget runs out:
        # that decides nothing, so it must not be reported as "no escape".
        print("  round %d: UNDECIDED -- the solver ran out of its conflict "
              "budget before finding an escape or proving that none exists."
              % rounds, flush=True)
        break
    if not got:
        print("", flush=True)
        print("  round %d: NO ESCAPE SURVIVES.  Forcing these %d pairs apart "
              "BLOCKS the point." % (rounds, len(cover)), flush=True)
        for (a, b) in cover:
            d2 = (g.vertices[a] - g.vertices[b]).norm2()
            print("     %d,%d   d^2 = %s   (%.4f)" % (a, b, d2, float(d2)),
                  flush=True)
        print("  carrier + p would be 5-chromatic on %d vertices   [%.0fs]"
              % (n + 1, time.time() - t0), flush=True)
        break
    sets = []
    for col in got:
        mono = set()
        for i, a in enumerate(NB):
            for b in NB[i + 1:]:
                if col[a] == col[b] and (a, b) not in adj and (b, a) not in adj:
                    mono.add((a, b))
        sets.append(mono)
    if any(not m for m in sets):
        print("  round %d: an escape has no monochromatic non-adjacent pair -- "
              "no family of forced pairs can block this point" % rounds,
              flush=True)
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
    print("  round %d: %d new escapes, +%d pairs, cover now %d   [%.0fs]"
          % (rounds, len(got), added, len(cover), time.time() - t0), flush=True)
json.dump({"n": n, "neighbourhood": NB, "rounds": rounds,
           "escapes_seen": total_escapes,
           "cover": [[a, b, str((g.vertices[a] - g.vertices[b]).norm2())]
                     for a, b in cover]},
          open(ROOT + "/data/blocking_requirements.json", "w"))
print("  written data/blocking_requirements.json (%d pairs)" % len(cover),
      flush=True)
