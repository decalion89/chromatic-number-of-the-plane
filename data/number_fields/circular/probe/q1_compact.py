"""Compact Q1 tables (one line per symmetry orbit of components)."""
from snr import *
R_LIST = [F(3001, 10000), F(302, 1000), F(305, 1000), F(31, 100), F(32, 100), F(1, 3)]
from collections import defaultdict
for k in (1, 2, 3):
    N = 5 ** k
    print(f"\n#### N = {N}\n")
    print("| r | #comp | total area | orbit (size) | kappa | point of the component (kappa attained) | area | nearest type pt | max dist |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in R_LIST:
        L = lift_components(k, r)
        tot = sum(area(P) for P in L)
        orbits = defaultdict(list)
        for P in L:
            orbits[orbit_key(P, N)].append(clean(P))
        rows = []
        for key, Ps in orbits.items():
            P = Ps[0]
            kap, opt, how = component_kappa_shifted(P, k)
            name, TP, dmax, dmin = classify(P, k)
            rows.append((kap, len(Ps), opt, area(P), name, TP, dmax))
        rows.sort(key=lambda t: (-t[0], t[4], float(t[6])))
        first = True
        for kap, n, opt, ar, name, TP, dmax in rows:
            lab = "C" if kap == F(1, 2) else ("Q" if kap == F(1, 3) else "X")
            tp = "c*" if name == "c" else "q"
            pt = f"({opt[0]}, {opt[1]})" if lab == "X" else ("c*" if lab == "C" else "q")
            rr = f"{r}" if first else ""
            nc = f"{len(L)}" if first else ""
            ta = f"{float(tot):.6f}" if first else ""
            print(f"| {rr} | {nc} | {ta} | {lab} ({n}) | {kap} | {pt} | {ar if ar.denominator < 10**7 else float(ar)} ~ {float(ar):.3e} | {tp} | {math.sqrt(dmax):.5f} |")
            first = False
