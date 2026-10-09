"""lift enumeration at a fixed r up to level K, with genealogy.  For each level: number of components, which are main
(contain c* or a point of Q_N, tested exactly), and for the others ('X') their exact kappa (certified LP), the type of
their parent (C, Q or X) and the max distance to the nearest type point.  r0(k) = max kappa of the X components."""
from snr import *
import sys, time
from collections import Counter
r = F(sys.argv[1]); K = int(sys.argv[2])
lo, hi = r, 1 - r
comps = [([(lo, lo), (hi, lo), (hi, hi), (lo, hi)], "R")]
for kk in range(1, K + 1):
    t0 = time.time()
    Np = 5 ** (kk - 1); N = 5 ** kk
    fams = []
    for j in (kk, -kk):
        A, B = rho_pow(j); fams += [(A, B), (B, -A)]
    new = []
    for C, lab in comps:
        for a in range(5):
            for b in range(5):
                pieces = [[(p[0] + Np * a, p[1] + Np * b) for p in C]]
                for (fa, fb) in fams:
                    nxt = []
                    for P in pieces:
                        nxt.extend(q for _, q in strip_split(P, fa, fb, lo, hi))
                    pieces = nxt
                    if not pieces: break
                new.extend((P, lab) for P in pieces)
    cs, qs = type_points(kk)
    funcs = functionals(kk)
    out = []; cnt = Counter(); kap_other = Counter(); maxd = F(0); mind = None
    for P, plab in new:
        pt = centroid_pt(P)
        lab = "X"
        for name, TP in [("C", cs[0])] + [("Q", q) for q in qs]:
            if all(floor(fa * TP[0] + fb * TP[1] - lo) == floor(fa * pt[0] + fb * pt[1] - lo) for (_, _, (fa, fb)) in funcs):
                lab = name
        if lab == "X":
            kap = component_kappa_fast(P, kk)[0]
            kap_other[(kap, plab)] += 1
            name, TP, dmax, dmin = classify(P, kk)
            maxd = max(maxd, dmax)
            mind = dmin if mind is None else min(mind, dmin)
        cnt[lab] += 1
        out.append((P, lab))
    comps = out
    print(f"level {kk} N={N}: {len(comps)} comps {dict(cnt)}; X at distance in [{math.sqrt(mind) if mind else 0:.2f}, {math.sqrt(maxd):.2f}] from type points ({time.time()-t0:.1f}s)")
    for (kv, pl) in sorted(kap_other, key=lambda t: (-t[0], t[1])):
        print(f"     X kappa {kv} = {float(kv):.6f} parent {pl}: {kap_other[(kv, pl)]}")
    sys.stdout.flush()
