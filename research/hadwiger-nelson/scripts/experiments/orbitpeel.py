"""Peel by orbits, not by vertices: the tightest carrier that still has rings.

Two requirements pull against each other and this work has recorded both
separately without ever putting them together.

A cap needs a TIGHT carrier: chromatic number equal to the number of colours
asked about, or a free colour makes every ring showable.  Tightening means
peeling towards criticality.

A cap also needs RINGS: a centre with a full orbit of points around it, which
is what makes the palette bound bite on something.  And peeling towards
criticality destroys rings -- minimising G's disjunction witness from 1581
vertices to 63 kept only what was load-bearing for the disjunction, and ring
structure was not load-bearing for that, so it went.

The way out is to peel in units that cannot break a ring.  The group here is
the order-twelve dihedral group of the hexagonal rotation and the mirror, and
a ring is a union of its orbits.  Removing a whole orbit removes ring points
in complete symmetric sets or not at all; it never leaves a ring with a hole.
So: order the orbits, try to drop each in turn, keep the drop when what is
left still needs five colours, and stop when nothing more can go.

What comes out is the tightest 5-chromatic carrier reachable from G without
ever breaking its symmetry -- which is the carrier the cap test has been
missing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
t0 = time.time()
builders = {"Sa": build_Sa, "Y": build_Y,
            "G": lambda f: build_G(f, as_graph=False)}
P = list(builders[SEED](K))
n = len(P)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
print(f"{SEED}: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]",
      flush=True)

# Orbits under the order-twelve group about the origin.
rot = _rot60(K)
idx = {p: i for i, p in enumerate(P)}
orbit = [None] * n
no = 0
for i, p in enumerate(P):
    if orbit[i] is not None:
        continue
    members = set()
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            if q in idx:
                members.add(idx[q])
            q = rot(q)
    for j in members:
        if orbit[j] is None:
            orbit[j] = no
    no += 1
groups = defaultdict(list)
for i, o in enumerate(orbit):
    groups[o].append(i)
sizes = defaultdict(int)
for o, v in groups.items():
    sizes[len(v)] += 1
print(f"{no} orbits under the order-12 group; sizes {dict(sorted(sizes.items()))}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def chromatic_at_least(keep, kk):
    """True when the induced subgraph on `keep` needs more than kk colours."""
    ks = set(keep)
    loc = {v: i for i, v in enumerate(sorted(ks))}
    m = len(loc)
    cls = [[1 + i * kk + c for c in range(kk)] for i in range(m)]
    for a, c in E:
        if a in ks and c in ks:
            for col in range(kk):
                cls.append([-(1 + loc[a] * kk + col), -(1 + loc[c] * kk + col)])
    for col in range(1, kk):
        cls.append([-(1 + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    out = s.solve()
    s.delete()
    return not out


keep = set(range(n))
assert chromatic_at_least(keep, k - 1), f"{SEED} does not need {k} colours"
print(f"confirmed: {SEED} needs {k} colours  [{time.time()-t0:.0f}s]",
      flush=True)
order = sorted(groups, key=lambda o: -len(groups[o]))
dropped = 0
for o in order:
    cand = keep - set(groups[o])
    if not cand:
        continue
    if chromatic_at_least(cand, k - 1):
        keep = cand
        dropped += 1
        if dropped % 10 == 0:
            print(f"   dropped {dropped} orbits, {len(keep)} points left"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
kept_orbits = len({orbit[v] for v in keep})
sub = sorted(keep)
ke = [(a, c) for a, c in E if a in keep and c in keep]
print(f"\norbit-peeled {SEED}: {len(sub)} points, {len(ke)} edges, "
      f"{kept_orbits} orbits, still {k}-chromatic  [{time.time()-t0:.0f}s]",
      flush=True)
print(f"dropped {dropped} of {no} orbits", flush=True)

# What rings survive?  A ring is a same-distance class around a kept centre.
adj = defaultdict(set)
for a, c in ke:
    adj[a].add(c)
    adj[c].add(a)
deg = {v: len(adj[v]) for v in sub}
dm, D2 = b.dim, b.D * b.D
import numpy as np
best = sorted(sub, key=lambda v: -deg[v])[:12]
for ci in best[:6]:
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    g = defaultdict(int)
    for off in np.nonzero(ok)[0]:
        if int(off) in keep and off != ci:
            D = Fr(int(sq[off, 0]), D2)
            if 0 < D <= 4:
                g[D] += 1
    big = sorted(((c, str(D)) for D, c in g.items()), reverse=True)[:5]
    print(f"   centre {ci} (deg {deg[ci]}): rings {big}", flush=True)
import pickle
with open(f"/tmp/claude-0/-home-user-darwin-50/"
          f"aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
          f"peeled_{SEED}_{k}.pkl", "wb") as f:
    pickle.dump([P[i] for i in sub], f)
print("DONE", flush=True)
