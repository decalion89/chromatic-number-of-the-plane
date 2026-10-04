"""follow components with kappa >= r through levels (lift enumeration), print positions/values of non-main ones"""
from snr import *
import sys, time
r = F(sys.argv[1]); K = int(sys.argv[2])
lo, hi = r, 1 - r
comps = [[(lo, lo), (hi, lo), (hi, hi), (lo, hi)]]
for kk in range(1, K + 1):
    t0 = time.time()
    Np = 5 ** (kk - 1); N = 5 ** kk
    fams = []
    for j in (kk, -kk):
        A, B = rho_pow(j); fams += [(A, B), (B, -A)]
    new = []
    for C in comps:
        for a in range(5):
            for b in range(5):
                pieces = [[(p[0] + Np * a, p[1] + Np * b) for p in C]]
                for (fa, fb) in fams:
                    nxt = []
                    for P in pieces:
                        nxt.extend(q for _, q in strip_split(P, fa, fb, lo, hi))
                    pieces = nxt
                    if not pieces: break
                new.extend(pieces)
    comps = new
    print(f"level {kk} N={N}: {len(comps)} components ({time.time()-t0:.1f}s)")
    for P in comps:
        name, TP, dmax, dmin = classify(P, kk)
        if dmax < 1: continue
        pt = P[0] if len(P) <= 2 else centroid_pt(P)
        print(f"   other: nverts={len(P)} pt=({pt[0]},{pt[1]}) ~({float(pt[0]):.4f},{float(pt[1]):.4f}) c/N~({float(pt[0]/N):.5f},{float(pt[1]/N):.5f}) near {name} dist {math.sqrt(dmax):.3f} kappa_pt={kappa_point(pt,kk)}")
