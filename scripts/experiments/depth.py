"""How fast does gluing shrink the interface, and how deep would five need?

At four colours the glue halves the space of interface patterns -- Sa carries
72 on eight high-degree vertices, and Sa glued to a rotated copy carries 38 --
and that shrinking is the mechanism: it continues until pairs have nowhere left
to disagree.  At five colours the same eight vertices carry 653 patterns in G
and 921 in the 807, with no saturation at all.

So the question is not whether gluing helps at five but at what RATE.  Glue
repeatedly, always at the highest-overlap glue available, and measure the same
eight vertices each time.  A constant factor per glue turns "it does not move"
into a number of levels; no factor at all turns it into a different kind of
statement.

The interface vertices are pinned by coordinate, not by index, so the count is
comparable across depths.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
from fractions import Fraction as Fr
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

SRC, K, DEPTH, M = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
if SRC == "Sa":
    F = Field((3, 11, 247)); pts = build_Sa(F)
elif SRC == "G":
    F = Field((3, 5, 7, 11)); pts = build_G(F, as_graph=False)
else:
    d = json.load(open(f"{HN_DIR}/data/{SRC}"))
    F = Field(tuple(d["field_generators"]))
    pts = [Point(F.element([Fr(a, b) for a, b in x]),
                 F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(pts)
TPTS = [g0.vertices[v] for v in sorted(range(g0.n), key=lambda v: -len(g0.adj[v]))[:8]]
ROTS = [("rot60", rotation_joining(Fr(1), F)),
        ("rot120", rotation_joining(Fr(1, 3), F)),
        ("rot180", rotation_joining(Fr(1, 4), F))]

def canon(p):
    seen, out = {}, []
    for c in p:
        if c not in seen: seen[c] = len(seen)
        out.append(seen[c])
    return tuple(out)

def count_patterns(pts, M):
    g = build_graph(pts); n = g.n
    pos = {p: i for i, p in enumerate(g.vertices)}
    T = [pos[p] for p in TPTS if p in pos]
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    for v in range(n):                 # at-most-one, so the readout is exact
        for c in range(K):
            for e in range(c + 1, K):
                cnf.append([-X(v, c), -X(v, e)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve(): s.delete(); return None, n, g, len(T)
    rng = random.Random(9)
    pos0 = set(l for l in s.get_model() if l > 0)
    col = [next(c for c in range(K) if X(v, c) in pos0) for v in range(n)]
    seen = {canon([col[v] for v in T])}
    for i in range(M):
        samp = rng.sample(range(n), 40)
        s.add_clause([-X(v, col[v]) for v in samp])
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): break
        p2 = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in p2) for v in range(n)]
        seen.add(canon([col[v] for v in T]))
    s.delete()
    return len(seen), n, g, len(T)

t0 = time.time()
cur = pts
for depth in range(DEPTH + 1):
    cnt, n, g, tsz = count_patterns(cur, M)
    m = sum(len(a) for a in g.adj) // 2
    print(f"  depth {depth}: n={n} m={m} deg={2.0*m/n:.2f}  interface {tsz} vertices"
          f"  patterns {cnt}   [{time.time()-t0:.0f}s]", flush=True)
    if cnt is None:
        print(f"  *** depth {depth} refuses {K} colours ***", flush=True); break
    if depth == DEPTH: break
    S = set(cur); best = None
    for name, r0 in ROTS:
        for c in cur[:600]:
            rot = r0.about(c)
            ov = sum(1 for p in cur if rot(p) in S)
            if ov < len(cur) and (best is None or ov > best[0]): best = (ov, rot)
    ov, rot = best
    seen2, nxt = set(), []
    for p in cur:
        for q in (p, rot(p)):
            if q not in seen2: seen2.add(q); nxt.append(q)
    print(f"    glue overlap {ov}/{len(cur)} -> {len(nxt)}", flush=True)
    cur = nxt
