"""Exact 3D configurations: a point is q / sqrt2 with q in Q^3 (rationals as strings in JSON); squared distances are
|q - q'|^2 / 2, rational. Graph, independent sets and congruence relations computed exactly."""
import json, sys
from fractions import Fraction as F
from collections import defaultdict
import numpy as np

def from_npy(path, den=100000):
    X = np.load(path) * np.sqrt(2); pts = []
    for p in X:
        q = [F(float(x)).limit_denominator(den) for x in p]
        assert max(abs(float(c) - x) for c, x in zip(q, p)) < 1e-10, ('not rational over sqrt2', p)
        pts.append(q)
    return pts
def save(path, pts): json.dump(dict(scale='1/sqrt2', points=[[str(c) for c in q] for q in pts]), open(path, 'w'))
def load(path):
    d = json.load(open(path)); assert d['scale'] == '1/sqrt2'
    return [[F(c) for c in q] for q in d['points']]
def d2(p, q): return sum((a - b) ** 2 for a, b in zip(p, q)) / 2
def floats(pts): return np.array([[float(c) for c in q] for q in pts]) / np.sqrt(2)
def graph(pts):
    n = len(pts); D2 = [[d2(pts[i], pts[j]) for j in range(n)] for i in range(n)]
    assert all(D2[i][j] != 0 for i in range(n) for j in range(n) if i != j), 'repeated point'
    return D2, [(i, j) for i in range(n) for j in range(i + 1, n) if D2[i][j] == 1]
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
def congruences(n, D2, atoms, minsize=3):
    """[(I, J)]: each non-representative member J of a congruence class of independent sets (>= minsize points)
    against the representative I, listed in matching order (all squared distances equal, exactly)"""
    def sig(S): return tuple(sorted(D2[a][b] for a in S for b in S if a < b))
    def bij(S, T):
        k = len(S); used = [False] * k; img = [None] * k
        def bt(p):
            if p == k: return True
            for c in range(k):
                if used[c]: continue
                if all(D2[S[q]][S[p]] == D2[T[img[q]]][T[c]] for q in range(p)):
                    used[c] = True; img[p] = c
                    if bt(p + 1): return True
                    used[c] = False
            return False
        return [T[c] for c in img] if bt(0) else None
    groups = defaultdict(list)
    for m in atoms:
        S = tuple(i for i in range(n) if m >> i & 1)
        if len(S) >= minsize: groups[(len(S), sig(S))].append(S)
    pairs = []; ncl = 0
    for mem in groups.values():
        reps = []
        for S in mem:
            for rep, lst in reps:
                b = bij(rep, S)
                if b is not None: lst.append(b); break
            else: reps.append((S, []))
        for rep, lst in reps:
            ncl += 1
            for T in lst: pairs.append((list(rep), T))
    return pairs, ncl
if __name__ == '__main__':
    save(sys.argv[2], from_npy(sys.argv[1]))
