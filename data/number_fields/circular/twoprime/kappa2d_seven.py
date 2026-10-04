"""As kappa2d.py, but labels the 7-adic components too: a component is '7' if one of the points N(a+bi)/7
(N = 5^k 13^m, a, b mod 7) has the same strip indices.  Prints the exact kappa of every component that is not
c, q or 7 (certified LP), grouped.  usage: python3 kappa2d_seven.py r schedule"""
import sys, time
from collections import Counter
from probe2d import *
from kappa2d import comp_kappa
r = Fr(sys.argv[1]); sched = sys.argv[2].split(',')
P = Probe(r)
t0 = time.time()
for step in sched:
    P.lift5() if step == '5' else P.lift13()
    print(f"  ({P.k},{P.m}): {len(P.comps)} components ({time.time()-t0:.0f}s)", flush=True)
N = P.N
fams = P.all_families()
lo = P.r
sevens = [(Fr(N * a, 7), Fr(N * b, 7)) for a in range(7) for b in range(7) if (a, b) != (0, 0)]
def idx(p): return tuple(floor(fa * p[0] + fb * p[1] - lo) for (fa, fb) in fams)
seven_idx = {}
for p in sevens:
    if all(lo <= (fa * p[0] + fb * p[1]) - floor(fa * p[0] + fb * p[1]) <= 1 - lo for (fa, fb) in fams):
        seven_idx[idx(p)] = p
lab = Counter(); kap = Counter(); ex = {}
for C in P.comps:
    pt = centroid(C)
    L = P.label(pt)
    if L == 'X' and idx(pt) in seven_idx:
        L = '7'
    lab[L] += 1
    if L == 'X':
        t, opt = comp_kappa(P, C, P.k, P.m)
        kap[t] += 1; ex.setdefault(t, opt)
print(f"r = {r}, (k,m)=({P.k},{P.m}), N={N}: {len(P.comps)} comps {dict(lab)}; 7-torsion points in S: {len(seven_idx)}")
for t in sorted(kap, reverse=True)[:30]:
    print(f"   other: kappa = {t} = {float(t):.7f}: {kap[t]} components; optimum {ex[t]}")
print("max kappa of the other components:", max(kap) if kap else None)
