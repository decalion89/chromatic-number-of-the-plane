"""Two forced pairs sharing a pivot give a spindle with a free angle.

The usual spindle needs |a-b| to satisfy that 4|a-b|^2 - 1 is a square, because
the rotation has to carry b to distance exactly 1 from ITSELF.  But forcing is
transitive, so the forced-equal pairs fall into classes, and a class with three
members a, b, c gives something better: any isometry rho fixing a with
|rho(b) - c| = 1 is enough, since then

    c(rho(b)) = c(rho(a)) = c(a) = c(c)   and   rho(b) ~ c.

The rotation no longer has to close a circle on itself; it only has to bring b
to somewhere on the unit circle about c.  That exists whenever

    | |a-b| - |a-c| |  <=  1  <=  |a-b| + |a-c|

and its cosine is fixed by the triangle -- one value, computable, and the only
question is whether it lies in the field.  A whole extra degree of freedom over
the classical spindle, available as soon as a forced class has three members.

This computes the classes at each level of the densification chain and reports
any class of size three or more together with the angles it offers.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

LEVELS = int(sys.argv[1]) if len(sys.argv) > 1 else 4
F = Field((3, 11, 247))
K = 4
REP = {1, 3, 11, 33, 247, 741, 2717, 8151}

def sqfree(r):
    d = 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r

def pairs_of(pts):
    g = build_graph(pts); n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve(): s.delete(); return None, g
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
    r = random.Random(5)
    for _ in range(72):
        if len(cols) >= 24: break
        s.set_phases([(1 if r.random() < 0.25 else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        p2 = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
    s.delete()
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(c[v] for c in cols)].append(v)
    out = []
    for vs in buck.values():
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s2 = Solver(name="m22", bootstrap_with=cnf)
                d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
                if not d: out.append((a, b))
    return out, g

def classes(pairs):
    par = {}
    def find(x):
        par.setdefault(x, x)
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb
    out = defaultdict(list)
    for x in par: out[find(x)].append(x)
    return [sorted(v) for v in out.values() if len(v) > 1]

pts = build_Sa(F)
t0 = time.time()
for lvl in range(LEVELS):
    pr, g = pairs_of(pts)
    cls = classes(pr) if pr else []
    sizes = Counter(len(c) for c in cls)
    print(f"\nlevel {lvl}: n={g.n}  forced pairs {len(pr) if pr else 0}  "
          f"classes {dict(sorted(sizes.items()))}   [{time.time()-t0:.0f}s]", flush=True)
    for c in cls:
        if len(c) < 3: continue
        print(f"  class of {len(c)}: {c[:8]}", flush=True)
        found = 0
        for ai in c:
            for bi in c:
                for ci in c:
                    if len({ai, bi, ci}) < 3: continue
                    ab = (pts[ai] - pts[bi]).norm2()
                    ac = (pts[ai] - pts[ci]).norm2()
                    if not (ab.is_rational() and ac.is_rational()): continue
                    r1, r2 = Fr(ab.c[0]), Fr(ac.c[0])
                    # a rotation about a takes b to the unit circle about c iff
                    # the circles of radius |a-b| about a and 1 about c meet
                    lo = r1 + r2 - 2 * 1  # compare (|ab|-|ac|)^2 <= 1 numerically
                    import math
                    d1, d2 = math.sqrt(float(r1)), math.sqrt(float(r2))
                    if abs(d1 - d2) > 1 or d1 + d2 < 1: continue
                    # cos of the angle between (b-a) and (c-a) after rotation:
                    # |rho(b) - c|^2 = r1 + r2 - 2 sqrt(r1 r2) cos t = 1
                    num = r1 + r2 - 1
                    cos2 = Fr(num * num, 1) / (4 * r1 * r2) if r1 * r2 else None
                    if cos2 is None or cos2 > 1: continue
                    s2 = 1 - cos2
                    rad = 1 if s2 == 0 else sqfree(s2.numerator * s2.denominator)
                    if rad in REP:
                        found += 1
                        if found <= 4:
                            print(f"      a={ai} b={bi} c={ci}: |ab|^2={r1} "
                                  f"|ac|^2={r2}  cos^2={cos2}  radical {rad}"
                                  f"   <<<< FREE ANGLE", flush=True)
        print(f"      {found} usable triangles in this class", flush=True)
    if lvl + 1 >= LEVELS: break
    S = set(pts); best = None
    for name, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3))):
        r0 = rotation_joining(d2, F)
        for cpt in pts:
            rot = r0.about(cpt)
            ov = sum(1 for p in pts if rot(p) in S)
            if ov < len(pts) and (best is None or ov > best[0]): best = (ov, rot)
    ov, rot = best
    seen, nxt = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen: seen.add(q); nxt.append(q)
    print(f"  glue overlap {ov} -> {len(nxt)}", flush=True)
    pts = nxt
