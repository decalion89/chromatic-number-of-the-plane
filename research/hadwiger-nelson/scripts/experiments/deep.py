"""Densify at four colours as far as it will go, and watch the five-colour slack.

The measured trajectory that motivates it:

    level 0  Sa   397 pts  deg  9.94  free@4 0.00%   0 forced pairs
    level 1       570      deg 11.36  free@4 0.00%   8
    level 2       679      deg 12.05  free@4 0.00%  27

and the spindle of level 1 is a 1139-point 5-chromatic graph with free@5
11.33%.  A random graph of mean degree d has roughly k(1-1/k)^d free vertices,
which at k=4, d=10 would be 22%; Sa is at 0.00%, so its structure is worth
about a factor of three in degree.  Applying the same factor at k=5 puts the
threshold near degree 14 to 16.  The chain is at 12.05 and climbing while the
size growth decays, so the question is simply whether it gets there.

Two selection rules are run against each other: the glue with the largest
overlap, and the glue that actually yields the largest mean degree (evaluated
on the top candidates, since that needs the edges).
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

LEVELS = int(sys.argv[1]); RULE = sys.argv[2]; NEVAL = int(sys.argv[3])
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

def forced(g, pts2, cols, cnf, X):
    buck = defaultdict(list)
    for v in range(g.n):
        buck[tuple(col[v] for col in cols)].append(v)
    out = []
    for vs in buck.values():
        if len(vs) < 2: continue
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s = Solver(name="m22", bootstrap_with=cnf)
                d = s.solve(assumptions=[X(a, 0), X(b, 1)]); s.delete()
                if not d:
                    dd = (pts2[a] - pts2[b]).norm2()
                    out.append((a, b, dd,
                                radical(Fr(dd.c[0])) if dd.is_rational() else None))
    return out

def apply(pts, f):
    seen, out = set(), []
    for p in pts:
        for q in (p, f(p)):
            if q not in seen: seen.add(q); out.append(q)
    return out

def candidates(pts, howmany):
    S = set(pts); n0 = len(pts); out = []
    diff = Counter()
    for p in pts:
        for q in pts:
            if p is not q: diff[p - q] += 1
    for t, mult in diff.most_common(2 * howmany):
        out.append((mult, "translate", (lambda p, t=t: p + t)))
    for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4))):
        r0 = rotation_joining(d2, F)
        for c in pts:
            rot = r0.about(c)
            ov = sum(1 for p in pts if rot(p) in S)
            if ov < n0: out.append((ov, name, rot))
    out.sort(key=lambda x: -x[0])
    return out

pts = build_Sa(F)
t0 = time.time()
hist = []
for lvl in range(LEVELS):
    g = build_graph(pts)
    m = sum(len(a) for a in g.adj) // 2
    cols4, cnf4, X4 = colourings(g, 4, 12)
    sl4 = None if not cols4 else slack(g, 4, cols4)
    print(f"\n=== level {lvl}: n={g.n} m={m} deg={2.0*m/g.n:.2f} "
          f"4-col={bool(cols4)} free@4={'--' if sl4 is None else f'{sl4:.2f}%'}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not cols4:
        print("  *** the carrier refuses four colours ***", flush=True); break
    fp = forced(g, pts, cols4, cnf4, X4)
    good = [f for f in fp if f[3] in REP]
    print(f"  forced pairs at 4: {len(fp)} ({len(good)} spindleable); "
          f"distances {dict(Counter(str(f[2]) for f in fp).most_common(3))}", flush=True)
    if good:
        a, b, d2, rad = good[0]
        Zp = apply(pts, rotation_joining(Fr(d2.c[0]), F).about(pts[a]))
        gz = build_graph(Zp); mz = sum(len(x) for x in gz.adj) // 2
        c4, _, _ = colourings(gz, 4, 1)
        c5, cnf5, X5 = colourings(gz, 5, 12)
        sl5 = None if not c5 else slack(gz, 5, c5)
        print(f"  -> spindle: n={gz.n} m={mz} deg={2.0*mz/gz.n:.2f} 4-col={bool(c4)} "
              f"5-col={bool(c5)} free@5={'--' if sl5 is None else f'{sl5:.2f}%'}",
              flush=True)
        hist.append((lvl, g.n, 2.0*m/g.n, len(fp), gz.n, 2.0*mz/gz.n, sl5))
        if not c5:
            print("  *** SIX COLOURS ***", flush=True)
            pickle.dump((RULE, lvl), open(f"/tmp/hn/SIX_{RULE}_{lvl}.pkl", "wb")); break
        if not c4:
            f5 = forced(gz, Zp, c5, cnf5, X5)
            print(f"     forced pairs at FIVE: {len(f5)}"
                  + ("  <<<<<<<<<<<<<<" if f5 else ""), flush=True)
            for x in f5[:8]:
                print(f"        v{x[0]} v{x[1]} d^2={x[2]} radical {x[3]}", flush=True)
            if f5:
                pickle.dump((RULE, lvl, [(x[0], x[1], str(x[2]), x[3]) for x in f5]),
                            open(f"/tmp/hn/FIVE_{RULE}_{lvl}.pkl", "wb"))
    cand = candidates(pts, 60)
    if RULE == "overlap":
        ov, name, f = cand[0]
        nxt = apply(pts, f)
        print(f"  glue: {name} overlap {ov} -> {len(nxt)}", flush=True)
    else:
        bestd = None
        for ov, name, f in cand[:NEVAL]:
            trial = apply(pts, f)
            gt = build_graph(trial)
            mt = sum(len(x) for x in gt.adj) // 2
            d = 2.0 * mt / gt.n
            if bestd is None or d > bestd[0]: bestd = (d, ov, name, trial)
        d, ov, name, nxt = bestd
        print(f"  glue: {name} overlap {ov} deg {d:.2f} -> {len(nxt)}", flush=True)
    pts = nxt
print("\nlevel  n_4    deg_4   forced  n_5    deg_5   free@5", flush=True)
for h in hist:
    print(f"{h[0]:>5} {h[1]:>6} {h[2]:>7.2f} {h[3]:>7} {h[4]:>6} {h[5]:>7.2f} "
          f"{'--' if h[6] is None else f'{h[6]:.2f}%'}", flush=True)
