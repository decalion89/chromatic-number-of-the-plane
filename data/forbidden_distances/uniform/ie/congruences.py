"""Congruence constraints for X23 generated from scratch: the independent sets of the exact unit-distance graph are
split into congruence classes (S ~ S' iff some bijection preserves all squared distances, which are compared exactly
in Q(sqrt3, sqrt11)); each class S_1, ..., S_m gives the m - 1 ordered pairs (S_1, S_j), with S_j listed in the order
of the bijection. Sets of size 1 are skipped (they repeat P(i in S) = delta). Writes cong23.json: [[I, J], ...] with
1-based indices. usage: python3 gen_cong.py [maxsize]"""
import sys, os, json, time
from collections import defaultdict
from kfield import K, ONE, dist2
from x23 import P
maxsize = int(sys.argv[1]) if len(sys.argv) > 1 else 99
n = len(P); pts = [(K(re), K(im)) for re, im in P]
d2 = [[None] * n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j: d2[i][j] = dist2(pts[i], pts[j])
nb = [sum(1 << j for j in range(n) if j != i and d2[i][j] == ONE) for i in range(n)]
ind = []
def rec(start, cur, mask, lst):
    if len(lst) >= 2: ind.append(tuple(lst))
    for v in range(start, n):
        if not (mask >> v) & 1 and len(lst) < maxsize: rec(v + 1, cur | (1 << v), mask | nb[v], lst + [v])
rec(0, 0, 0, [])
def sig(S):
    return tuple(sorted(hash(d2[a][b]) for a in S for b in S if a < b))
def bijection(S, T):
    """an ordering of T matching S's order with all squared distances equal, or None"""
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
t0 = time.time()
groups = defaultdict(list)
for S in ind: groups[(len(S), sig(S))].append(S)
pairs = []; nclass = 0
for key, mem in groups.items():
    reps = []                                   # [(representative, [members in matched order])]
    for S in mem:
        for rep, lst in reps:
            b = bijection(rep, S)
            if b is not None: lst.append(b); break
        else: reps.append((S, []))
    for rep, lst in reps:
        nclass += 1
        for T in lst: pairs.append(([a + 1 for a in rep], [b + 1 for b in T]))
print(f"{len(ind)} independent sets of size >= 2, {nclass} congruence classes, {len(pairs)} congruence constraints, "
      f"{time.time() - t0:.0f}s")
json.dump(pairs, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cong23.json'), 'w'))
