"""trace parents of components level by level (lift enumeration), for a given r"""
from snr import *
import sys
r = F(sys.argv[1]); K = int(sys.argv[2])
lo, hi = r, 1 - r
comps = [([(lo, lo), (hi, lo), (hi, hi), (lo, hi)], "root", None)]
for kk in range(1, K + 1):
    Np = 5 ** (kk - 1); N = 5 ** kk
    fams = []
    for j in (kk, -kk):
        A, B = rho_pow(j); fams += [(A, B), (B, -A)]
    new = []
    for idx, (C, lab, _) in enumerate(comps):
        for a in range(5):
            for b in range(5):
                pieces = [[(p[0] + Np * a, p[1] + Np * b) for p in C]]
                for (fa, fb) in fams:
                    nxt = []
                    for P in pieces:
                        nxt.extend(q for _, q in strip_split(P, fa, fb, lo, hi))
                    pieces = nxt
                    if not pieces: break
                for P in pieces:
                    new.append((P, None, (idx, lab, (a, b))))
    # label
    out = []
    for P, _, par in new:
        kap, opt, how = component_kappa(P, kk)
        name, TP, dmax, dmin = classify(P, kk)
        main = (name == "c" and kap == F(1, 2)) or (name.startswith("q") and kap == F(1, 3) and dmax < F(1, 100))
        lab = ("C" if name == "c" else "Q") if main and (name=="c" or kap==F(1,3)) else "X"
        out.append((P, lab, par, kap, name, dmax, opt))
    print(f"level {kk} (N={N}): {len(out)} components")
    from collections import Counter
    cnt = Counter((lab, par[1], str(kap)) for P, lab, par, kap, name, dmax, opt in out)
    for key in sorted(cnt): print("   ", key, cnt[key])
    comps = [(P, lab, par) for P, lab, par, kap, name, dmax, opt in out]
