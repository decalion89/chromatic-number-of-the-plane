"""How many directions can a critical graph carry?

Blocking is a statement about how the edge directions sit mod 5, and covering
PG(r-1,5) needs at least six hyperplanes.  The Moser spindle has 14 projective
directions in a rank-4 module and fails.  So the question "can a critical
unit-distance graph block" turns first on a measurable one:

    how rich can the direction set of a 4-critical subgraph be?

Greedy deletion returns a MINIMAL subgraph, and which one depends entirely on
the order, so a single run says nothing.  This samples many random orders on
Sa -- 397 points, 4-chromatic, rank-16 multiquadratic module, the richest
4-chromatic graph here -- and records the distribution of core size, direction
count and module rank.  If every core is spindle-sized whatever the order,
critical graphs are direction-poor and blocking one is hopeless for a reason,
not for want of searching.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, random, time
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.graph import build_graph
from hn.homcol import _coords, has_homomorphism
from pysat.solvers import Solver

g = build_graph(build_Sa())
E = list(g.edges())
n = len(g.vertices)
adj = [[] for _ in range(n)]
for a, b in E:
    adj[a].append(b); adj[b].append(a)
print(f"Sa: {n} points, {len(E)} edges", flush=True)

vec = {}
for a, b in E:
    d = (g.vertices[b].x - g.vertices[a].x, g.vertices[b].y - g.vertices[a].y)
    vec[(a, b)] = _coords(d[0]) + _coords(d[1])


def three_col(keep):
    idx = {v: i for i, v in enumerate(keep)}
    cls = [[1 + i * 3 + c for c in range(3)] for i in range(len(idx))]
    for a, b in E:
        if a in idx and b in idx:
            for c in range(3):
                cls.append([-(1 + idx[a] * 3 + c), -(1 + idx[b] * 3 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def profile(keep):
    raw = [vec[(a, b)] for a, b in E if a in keep and b in keep]
    den = 1
    for v in raw:
        for q in v:
            den = den * Fraction(q).denominator // gcd(den,
                                                       Fraction(q).denominator)
    ints = {tuple(int(Fraction(q) * den) for q in v) for v in raw}
    cont = 0
    for v in ints:
        for x in v:
            cont = gcd(cont, abs(x))
    if cont > 1:
        ints = {tuple(x // cont for x in v) for v in ints}
    rows = [list(map(Fraction, v)) for v in ints]
    rank, piv = 0, 0
    cols = len(rows[0]) if rows else 0
    for c in range(cols):
        r = next((i for i in range(rank, len(rows)) if rows[i][c]), None)
        if r is None:
            continue
        rows[rank], rows[r] = rows[r], rows[rank]
        f = rows[rank][c]
        rows[rank] = [x / f for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][c]:
                k = rows[i][c]
                rows[i] = [x - k * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return sorted(ints), rank


best, seen, t0 = None, {}, time.time()
rng = random.Random(20260920)
for trial in range(1, 400):
    order = list(range(n))
    rng.shuffle(order)
    keep = set(range(n))
    for v in order:
        keep.discard(v)
        if three_col(keep):
            keep.add(v)
    dirs, rank = profile(keep)
    key = (len(keep), len(dirs), rank)
    seen[key] = seen.get(key, 0) + 1
    if best is None or len(dirs) > best[1]:
        best = (len(keep), len(dirs), rank)
        phi, _ = has_homomorphism(dirs, 5)
        print(f"  trial {trial}: core {len(keep)} points, {len(dirs)} "
              f"directions, rank {rank}, "
              + ("BLOCKS" if phi is None else "coset colouring exists")
              + f"  [{time.time()-t0:.0f}s]", flush=True)
    if trial % 50 == 0:
        print(f"  {trial} trials, profiles seen: "
              + ", ".join(f"{k[0]}pts/{k[1]}dir/rk{k[2]}x{v}"
                          for k, v in sorted(seen.items())), flush=True)
