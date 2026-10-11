"""Exact configurations in Q(sqrt3, sqrt11): identification of floating-point points, the exact unit-distance graph,
the independent sets and the congruence relations between them (as in congruences.py, for any point set).
A point is ((a, b, c, d), (a', b', c', d')): real and imaginary parts a + b sqrt3 + c sqrt11 + d sqrt33."""
import json, time, os, sys
from fractions import Fraction as F
from collections import defaultdict
import mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'uniform', 'ie'))
from kfield import K, ONE, dist2

mp.mp.dps = 50
S3, S11, S33 = mp.sqrt(3), mp.sqrt(11), mp.sqrt(33)

def kval(k):
    a, b, c, d = k.v
    return mp.mpf(a.numerator) / a.denominator + mp.mpf(b.numerator) / b.denominator * S3 + \
        mp.mpf(c.numerator) / c.denominator * S11 + mp.mpf(d.numerator) / d.denominator * S33

def identify(x, maxcoeff=10 ** 4):
    """the element of K closest to x (a float or mpf) found by PSLQ, or None"""
    x = mp.mpf(x)
    if abs(x) < mp.mpf(10) ** -30: return K((0, 0, 0, 0))
    rel = mp.pslq([x, 1, S3, S11, S33], maxcoeff=maxcoeff, maxsteps=10 ** 6)
    if rel is None or rel[0] == 0: return None
    return K(tuple(F(-r, rel[0]) for r in rel[1:]))

def refine(pts_exact, z, tol=1e-7):
    """exact K point for the complex float z, using two points of pts_exact at distance 1 from z to recompute it in
    high precision before identifying it; checked exactly against those two points"""
    hp = [(kval(re), kval(im)) for re, im in pts_exact]
    nb = [i for i, (x, y) in enumerate(hp) if abs(abs(complex(float(x), float(y)) - z) - 1) < tol]
    assert len(nb) >= 2, 'not at distance 1 from two points'
    i, j = nb[0], nb[1]
    pi, pj = mp.mpc(*hp[i]), mp.mpc(*hp[j])
    d = pj - pi; r = abs(d); h = mp.sqrt(1 - r * r / 4); m = (pi + pj) / 2
    c = min((m + s * 1j * d / r * h for s in (1, -1)), key=lambda w: abs(complex(w) - z))
    re, im = identify(c.real), identify(c.imag)
    assert re is not None and im is not None, ('not in Q(sqrt3, sqrt11)', z)
    p = (re, im)
    for k in nb: assert dist2(p, pts_exact[k]) == ONE, 'identified point is not at distance 1'
    return p

def graph(pts):
    n = len(pts)
    d2 = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j: d2[i][j] = dist2(pts[i], pts[j])
    assert all(d2[i][j] != K((0, 0, 0, 0)) for i in range(n) for j in range(n) if i != j), 'repeated point'
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if d2[i][j] == ONE]
    return d2, edges

def independent_sets(n, edges):
    nb = [0] * n
    for i, j in edges: nb[i] |= 1 << j; nb[j] |= 1 << i
    out = []
    def rec(start, cur, mask):
        out.append(cur)
        for v in range(start, n):
            if not (mask >> v) & 1: rec(v + 1, cur | (1 << v), mask | nb[v])
    rec(0, 0, 0)
    return out

def congruences(n, d2, atoms):
    """[(I, J)] (0-based ordered lists): for each congruence class of independent sets with >= 2 points, its
    representative against every other member, matched by a bijection that preserves all squared distances"""
    sets = [tuple(i for i in range(n) if m >> i & 1) for m in atoms if m & (m - 1)]
    def sig(S): return tuple(sorted(hash(d2[a][b]) for a in S for b in S if a < b))
    def bijection(S, T):
        k = len(S); used = [False] * k; img = [None] * k
        def bt(p):
            if p == k: return True
            for c in range(k):
                if used[c]: continue
                if all(d2[S[q]][S[p]] == d2[T[img[q]]][T[c]] for q in range(p)):
                    used[c] = True; img[p] = c
                    if bt(p + 1): return True
                    used[c] = False
            return False
        return [T[c] for c in img] if bt(0) else None
    groups = defaultdict(list)
    for S in sets: groups[(len(S), sig(S))].append(S)
    pairs = []; ncl = 0
    for mem in groups.values():
        reps = []
        for S in mem:
            for rep, lst in reps:
                b = bijection(rep, S)
                if b is not None: lst.append(b); break
            else: reps.append((S, []))
        for rep, lst in reps:
            ncl += 1
            for T in lst: pairs.append((list(rep), T))
    return pairs, ncl

def kstr(k): return [str(x) for x in k.v]

def save(path, pts, extra=None):
    d = dict(points=[[kstr(re), kstr(im)] for re, im in pts])
    if extra: d.update(extra)
    json.dump(d, open(path, 'w'))

def load(path):
    d = json.load(open(path))
    return [(K(tuple(F(x) for x in re)), K(tuple(F(x) for x in im))) for re, im in d['points']], d
