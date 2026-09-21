"""Alternate the two minimisations until neither moves.

Two reductions have run separately and each helps the other.  Dropping pairs
frees vertices, since a vertex only matters through the pairs and edges that
touch it; dropping vertices frees pairs, since a pair with a removed endpoint
goes with it.  Run singly they stop early.

So alternate: minimise the forbidden pairs from several random orders and keep
the smallest result, then peel vertices, then go round again, until a full
round changes nothing.  On sixty-odd vertices each solve is milliseconds, so
many random orders are affordable -- which matters, because greedy removal
gives a MINIMAL set and the order decides which one.  The distance-descending
order gave 181 pairs and a random order gave 159, so the order is worth
sampling rather than trusting.

The output is meant to be small enough to print: a list of points with exact
coordinates, the unit-distance edges among them, and the forbidden pairs, such
that no proper 5-colouring avoids all of the latter.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
KEEP = [Fr(15, 16), Fr(16), Fr(17, 2), Fr(9), Fr(7)]
P = build_G(K, as_graph=False)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
n = len(P)
E = sorted(set((min(a, b), max(a, b))
               for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
byd = defaultdict(list)
for i in range(n - 1):
    d = rows[i + 1:] - rows[i]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    for off in np.nonzero(rat)[0]:
        j = i + 1 + int(off)
        if (i, j) not in Eset:
            v = int(sq[off, 0])
            if v:
                byd[Fr(v, D2)].append((i, j))

pairs = [(D, p) for D in KEEP for p in byd[D]]
verts = set()
for _, (a, b) in pairs:
    verts.add(a)
    verts.add(b)


def colours(vs, prs):
    vs = set(vs)
    ren = {v: i for i, v in enumerate(sorted(vs))}
    m = len(ren)
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, b in E:
        if a in vs and b in vs:
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    for _, (a, b) in prs:
        if a in vs and b in vs:
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


def min_pairs(vs, prs, tries=12):
    best = list(prs)
    for t in range(tries):
        cur = list(prs)
        rng = random.Random(500 + t)
        while True:
            snap = (sorted(cur, key=lambda x: -float(x[0])) if t == 0
                    else rng.sample(cur, len(cur)))
            dropped = 0
            for target in snap:
                if target not in cur:
                    continue
                trial = [x for x in cur if x is not target]
                if colours(vs, trial):
                    continue
                cur = trial
                dropped += 1
            if dropped == 0:
                break
        if len(cur) < len(best):
            best = cur
    return best


def min_verts(vs, prs):
    cur = set(vs)
    while True:
        dropped = 0
        for v in sorted(cur):
            if v not in cur:
                continue
            trial = cur - {v}
            if not colours(trial, prs):
                cur = trial
                dropped += 1
        if dropped == 0:
            return cur


print(f"start: {len(verts)} vertices, {len(pairs)} pairs"
      f"  [{time.time()-t0:.0f}s]", flush=True)
assert not colours(verts, pairs)
for rnd in range(1, 9):
    before = (len(verts), len(pairs))
    pairs = min_pairs(verts, pairs)
    pairs = [x for x in pairs if x[1][0] in verts and x[1][1] in verts]
    verts = min_verts(verts, pairs)
    pairs = [x for x in pairs if x[1][0] in verts and x[1][1] in verts]
    c = Counter(str(D) for D, _ in pairs)
    print(f"  round {rnd}: {len(verts)} vertices, {len(pairs)} pairs "
          f"{dict(sorted(c.items()))}  [{time.time()-t0:.0f}s]", flush=True)
    if (len(verts), len(pairs)) == before:
        break
ok = colours(verts, pairs)
print(f"\nfinal: {len(verts)} vertices, {len(pairs)} pairs, colours {ok}",
      flush=True)
deg = Counter()
for a, b in E:
    if a in verts and b in verts:
        deg[a] += 1
        deg[b] += 1
print(f"   unit edges among them: {sum(deg.values())//2}", flush=True)
print(f"   classes: {dict(sorted(Counter(str(D) for D, _ in pairs).items()))}",
      flush=True)
print("DONE", flush=True)
