"""Densify at four colours; at every level, spindle and try both routes to six.

The trajectory that motivates it:

    level 0  Sa       397 pts  deg  9.94  free@4 0.00%   0 forced pairs
    level 1           570      deg 11.36  free@4 0.00%   8 forced pairs
    level 2           679      deg 12.05  free@4 0.00%  27 forced pairs

Size growth decays (x1.44, then x1.19) while mean degree and the number of
forced pairs climb and the carrier stays perfectly rigid.  So the glue is a
densifier that costs less each time it is applied.

At each level the carrier is spindled at one of its forced pairs, which gives a
5-chromatic graph, and that graph is tried two ways:

  * a forced-equal pair at five, which one rotation would close (the two-copy
    route, de Grey's own);
  * three copies about a centre with a point at squared distance 1/3, whose
    images form a unit triangle -- that route needs only a palette of two
    rather than one, and four copies do not exist in the plane.

The glue alphabet is every translation by a frequent difference together with
the rotations about a vertex, ranked by overlap, because overlap is what
defeats the sigma-argument: with the copies disjoint, colouring the second by
sigma . c for a fixed-point-free sigma satisfies every cross edge at once.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

LEVELS = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NTRIPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 6
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
    if not s.solve():
        s.delete(); return [], cnf, X
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(g.n)]]
    rng = random.Random(seed); tries = 0
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
    return 100.0 * min(sum(1 for v in range(g.n)
                           if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
                       for col in cols) / g.n

def forced(g, pts2, cols, cnf, X, K):
    buck = defaultdict(list)
    for v in range(g.n):
        buck[tuple(col[v] for col in cols)].append(v)
    out = []
    for vs in buck.values():
        if len(vs) < 2: continue
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s = Solver(name="m22", bootstrap_with=cnf)
                differ = s.solve(assumptions=[X(a, 0), X(b, 1)]); s.delete()
                if not differ:
                    dd = (pts2[a] - pts2[b]).norm2()
                    out.append((a, b, dd,
                                radical(Fr(dd.c[0])) if dd.is_rational() else None))
    return out

def best_glues(pts, howmany):
    S = set(pts); n0 = len(pts); out = []
    diff = Counter()
    for p in pts:
        for q in pts:
            if p is not q: diff[p - q] += 1
    for t, mult in diff.most_common(4 * howmany):
        out.append((mult, "translate", (lambda p, t=t: p + t)))
    for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4))):
        r0 = rotation_joining(d2, F)
        for c in pts:
            rot = r0.about(c)
            ov = sum(1 for p in pts if rot(p) in S)
            if ov < n0: out.append((ov, name, rot))
    out.sort(key=lambda x: -x[0])
    return out

r120 = rotation_joining(Fr(1, 3), F)

def triple_scan(pts, K, howmany, tag):
    S = set(pts); n0 = len(pts)
    ranked = []
    for si, s in enumerate(pts):
        if not any((p - s).norm2() == Fr(1, 3) for p in pts): continue
        a = r120.about(s)
        seen = set()
        for p in pts:
            q1 = a(p); seen.add(p); seen.add(q1); seen.add(a(q1))
        if len(seen) == n0: continue
        ranked.append((3 * n0 - len(seen), si, a, len(seen)))
    ranked.sort(key=lambda x: -x[0])
    for ov, si, a, nn in ranked[:howmany]:
        seen, U = set(), []
        for p in pts:
            q1 = a(p); q2 = a(q1)
            for z in (p, q1, q2):
                if z not in seen: seen.add(z); U.append(z)
        g = build_graph(U)
        cnf, X = cnf_of(g, K)
        s = Solver(name="m22", bootstrap_with=cnf); ok = s.solve(); s.delete()
        print(f"     triple about v{si}: n={g.n} overlap={ov} "
              f"{K}-colourable={ok}" + ("" if ok else "   <<<<<<<<<<<<<<<<"),
              flush=True)
        if not ok:
            pickle.dump((tag, K, si), open(
                f"/tmp/hn/TRIPHIT_{tag}_{K}_{si}.pkl", "wb"))
            return True
    return False

pts = build_Sa(F)
t0 = time.time()
for lvl in range(LEVELS):
    g = build_graph(pts)
    m = sum(len(a) for a in g.adj) // 2
    cols4, cnf4, X4 = colourings(g, 4, 24)
    print(f"\n=== level {lvl}: n={g.n} m={m} deg={2.0*m/g.n:.2f} "
          f"4-colourable={bool(cols4)} free@4="
          f"{'--' if not cols4 else f'{slack(g,4,cols4):.2f}%'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not cols4:
        print("  *** the carrier itself refuses four colours ***", flush=True); break
    fp = forced(g, pts, cols4, cnf4, X4, 4)
    good = [f for f in fp if f[3] in REP]
    dd = Counter((str(f[2]), f[3]) for f in fp)
    print(f"  forced pairs at 4: {len(fp)}, spindleable {len(good)}; "
          f"distances {dict(list(dd.items())[:4])}", flush=True)
    if triple_scan(pts, 4, NTRIPLE, f"L{lvl}"):
        print("  *** the triple refuses four colours ***", flush=True)
    if good:
        a, b, d2, rad = good[0]
        rot = rotation_joining(Fr(d2.c[0]), F).about(pts[a])
        seen, Zp = set(), []
        for p in pts:
            for q in (p, rot(p)):
                if q not in seen: seen.add(q); Zp.append(q)
        gz = build_graph(Zp)
        mz = sum(len(x) for x in gz.adj) // 2
        c4, _, _ = colourings(gz, 4, 2)
        c5, cnf5, X5 = colourings(gz, 5, 24)
        print(f"  -> spindled at d^2={d2}: n={gz.n} m={mz} deg={2.0*mz/gz.n:.2f} "
              f"4-col={bool(c4)} 5-col={bool(c5)} free@5="
              f"{'--' if not c5 else f'{slack(gz,5,c5):.2f}%'}", flush=True)
        if not c4 and c5:
            f5 = forced(gz, Zp, c5, cnf5, X5, 5)
            print(f"     forced pairs at FIVE: {len(f5)}"
                  + ("   <<<<<<<<<<<<<<<<" if f5 else ""), flush=True)
            for x in f5[:6]:
                print(f"        v{x[0]} v{x[1]} d^2={x[2]} radical {x[3]}", flush=True)
            if f5:
                pickle.dump((lvl, [(x[0], x[1], str(x[2]), x[3]) for x in f5]),
                            open(f"/tmp/hn/FIVE_{lvl}.pkl", "wb"))
            triple_scan(Zp, 5, NTRIPLE, f"Z{lvl}")
        if not c5:
            print("  *** SIX COLOURS *** the spindled graph refuses five", flush=True)
            pickle.dump(lvl, open(f"/tmp/hn/SIX_{lvl}.pkl", "wb"))
            break
    glues = best_glues(pts, 40)
    ov, name, f = glues[0]
    seen, nxt = set(), []
    for p in pts:
        for q in (p, f(p)):
            if q not in seen: seen.add(q); nxt.append(q)
    print(f"  glue: {name}, overlap {ov}/{len(pts)} -> {len(nxt)} points", flush=True)
    pts = nxt
