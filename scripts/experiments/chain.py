"""Densify where tightness is reachable, and spindle once at the end.

At four colours Sa is perfectly tight -- every vertex's closed neighbourhood
already shows all four -- and gluing keeps it that way while RAISING the mean
degree: Sa is 9.94, and Sa glued to its 60-degree image about a vertex is
11.36, still at 0.00% free.  At five colours nothing is tight: G sits at 17%,
the new Z at 11.6%, and gluing Z again does not improve it.

So densify at four, where the carrier stays rigid, and spend the single
spindle at the end.  Each level:

    pick the glue with the largest overlap (keeps n down, degree up)
    check it is still 4-chromatic and measure the slack
    list the pairs forced to share a colour, and their spindle radicals
    spindle the best one and measure the five-colour slack of the result

If the five-colour slack falls as the four-colour carrier thickens, the
programme converges on a carrier that refuses five.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

LEVELS = int(sys.argv[1]) if len(sys.argv) > 1 else 4
F = Field((3, 11, 247))
REP = {1, 3, 11, 33, 247, 741, 2717, 8151}

def sqfree(r):
    d = 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r

def radical(q):
    if q <= Fr(1, 4): return None
    c = Fr(1) - Fr(1, 2) / q
    s2 = 1 - c * c
    if s2 == 0: return 1
    return sqfree(s2.numerator * s2.denominator)

ROTS = [("rot60", rotation_joining(Fr(1), F)),
        ("rot120", rotation_joining(Fr(1, 3), F)),
        ("moser", rotation_joining(Fr(3), F)),
        ("r5/9", rotation_joining(Fr(5, 9), F))]

def cnf_of(g, K):
    X = lambda v, c: 1 + v * K + c
    out = [[X(v, c) for c in range(K)] for v in range(g.n)]
    for u, v in g.edges():
        for c in range(K):
            out.append([-X(u, c), -X(v, c)])
    return out, X

def colourings(g, K, ncol, seed=3):
    cnf, X = cnf_of(g, K)
    s = Solver(name="m22", bootstrap_with=cnf)
    rng = random.Random(seed)
    cols, tries = [], 0
    while len(cols) < ncol and tries < ncol * 3:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(g.n) for c in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in pos) for v in range(g.n)])
    s.delete()
    return cols, cnf, X

def slack(g, K, cols):
    if not cols: return None
    best = 10 ** 9
    for col in cols:
        f = sum(1 for v in range(g.n)
                if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
        best = min(best, f)
    return 100.0 * best / g.n

def forced(g, pts2, K, cols, cnf, X):
    buck = defaultdict(list)
    for v in range(g.n):
        buck[tuple(col[v] for col in cols)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i+1:]]
    out = []
    for a, b in cand:
        s = Solver(name="m22", bootstrap_with=cnf)
        differ = s.solve(assumptions=[X(a, 0), X(b, 1)]); s.delete()
        if not differ:
            dd = (pts2[a] - pts2[b]).norm2()
            out.append((a, b, dd, radical(Fr(dd.c[0])) if dd.is_rational() else None))
    return out

pts = build_Sa(F)
t0 = time.time()
for lvl in range(LEVELS):
    g = build_graph(pts)
    m = sum(len(a) for a in g.adj) // 2
    cols4, cnf4, X4 = colourings(g, 4, 24)
    sl = slack(g, 4, cols4)
    print(f"\n=== level {lvl}: n={g.n} m={m} deg={2.0*m/g.n:.2f} "
          f"4-colourable={bool(cols4)} free@4={sl if sl is None else f'{sl:.2f}%'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not cols4:
        print("  *** the carrier itself refuses four colours ***", flush=True); break
    fp = forced(g, pts, 4, cols4, cnf4, X4)
    good = [f for f in fp if f[3] in REP]
    print(f"  forced-equal pairs at 4: {len(fp)};  spindleable in this field: {len(good)}",
          flush=True)
    seen_d = {}
    for f in fp:
        seen_d.setdefault((str(f[2]), f[3]), 0)
        seen_d[(str(f[2]), f[3])] += 1
    for (d, r), c in sorted(seen_d.items(), key=lambda kv: -kv[1])[:6]:
        print(f"      d^2={d:<28} radical {r:<6} x{c}", flush=True)

    # spindle the best pair and measure the five-colour slack of the result
    if good:
        a, b, dd, rad = good[0]
        rot = rotation_joining(Fr(dd.c[0]), F).about(pts[a])
        seen, Zp = set(), []
        for p in pts:
            for q in (p, rot(p)):
                if q not in seen: seen.add(q); Zp.append(q)
        gz = build_graph(Zp)
        mz = sum(len(x) for x in gz.adj) // 2
        c4, _, _ = colourings(gz, 4, 4)
        c5, cnf5, X5 = colourings(gz, 5, 24)
        s5 = slack(gz, 5, c5)
        print(f"  -> spindled at d^2={dd} : n={gz.n} m={mz} deg={2.0*mz/gz.n:.2f} "
              f"4-colourable={bool(c4)} free@5={s5 if s5 is None else f'{s5:.2f}%'}",
              flush=True)
        if not c4 and c5:
            f5 = forced(gz, Zp, 5, c5, cnf5, X5)
            g5 = [x for x in f5 if x[3] in REP]
            print(f"     chi=5 graph: forced-equal pairs at FIVE = {len(f5)}, "
                  f"spindleable = {len(g5)}", flush=True)
            for x in g5[:6]:
                print(f"        v{x[0]} v{x[1]} d^2={x[2]} radical {x[3]}  <<<<<<<<",
                      flush=True)
            pickle.dump((lvl, [(x[0], x[1], str(x[2]), x[3]) for x in f5]),
                        open(f"/tmp/hn/chain_{lvl}.pkl", "wb"))

    # next level: the glue with the largest overlap
    best = None
    for rname, rot0 in ROTS:
        for ci in range(len(pts)):
            rot = rot0.about(pts[ci])
            s = set(pts)
            ov = sum(1 for p in pts if rot(p) in s)
            if ov == len(pts): continue
            if best is None or ov > best[0]:
                best = (ov, rname, ci, rot)
    ov, rname, ci, rot = best
    seen, nxt = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen: seen.add(q); nxt.append(q)
    print(f"  glue chosen: {rname} about v{ci}, overlap {ov}/{len(pts)} "
          f"-> {len(nxt)} points", flush=True)
    pts = nxt
