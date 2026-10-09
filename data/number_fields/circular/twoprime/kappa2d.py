"""Exact kappa of every component of the two-prime probe S^r(k, m): the largest r' for which the component (same
strip indices) is still non-empty, an exactly certified LP.  Lists the distinct kappa values of the non-main (X)
components; the threshold of the window is the largest of them."""
import sys, time
from collections import Counter
from probe2d import *
from lpexact import lp_max

def comp_kappa(P, poly, k, m):
    pt = centroid(poly)
    fams = P.all_families(k, m)
    cons = []
    for (fa, fb) in fams:
        v = fa * pt[0] + fb * pt[1]
        n = floor(v)                    # pt lies in the strip [n + r, n + 1 - r] for this family
        base = v - n
        cons.append((fa, fb, base))          # t <= f(pt + d) - n
        cons.append((-fa, -fb, 1 - base))    # t <= n + 1 - f(pt + d)
    res = lp_max(cons)
    t, dx, dy, lam = res
    return t, (pt[0] + dx, pt[1] + dy)

if __name__ == "__main__":
    r = Fr(sys.argv[1]); sched = sys.argv[2].split(',')
    P = Probe(r)
    for step in sched:
        P.lift5() if step == '5' else P.lift13()
    lab = Counter()
    kap = Counter()
    pts = {}
    for C in P.comps:
        L = P.label(centroid(C))
        lab[L] += 1
        if L == 'X':
            t, opt = comp_kappa(P, C, P.k, P.m)
            kap[t] += 1
            pts.setdefault(t, opt)
    print(f"r = {r}, (k,m)=({P.k},{P.m}), N={P.N}: {len(P.comps)} comps {dict(lab)}")
    for t in sorted(kap, reverse=True):
        pt = pts[t]
        print(f"   X kappa = {t} = {float(t):.6f}: {kap[t]} components; an optimal point c = ({pt[0]}, {pt[1]})")
