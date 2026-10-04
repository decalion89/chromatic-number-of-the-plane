"""Question 1: S_N^(r) for N = 5, 25, 125 and r in R_LIST, exactly.  Two enumerations (lift, direct with symmetry;
for k <= 2 also direct without symmetry) are compared component by component (canonical vertex sets mod N).
For each orbit of components under the symmetry group (c -> ic, c -> conj c, translations by N Z[i]) prints:
kappa (exact), orbit size, a point (vertex average, exact), area (exact), the nearest type point and the maximal
distance of the component from it."""
from snr import *
import sys, time
from collections import defaultdict

R_LIST = [F(3001, 10000), F(302, 1000), F(305, 1000), F(31, 100), F(32, 100), F(1, 3)]
KS = [int(a) for a in sys.argv[1:]] or [1, 2, 3]


def keyset(polys, N):
    out = []
    for P in polys:
        Q = clean(P)
        c = centroid_pt(Q)
        shift = (N * floor(c[0] / N), N * floor(c[1] / N))
        out.append(tuple(sorted((v[0] - shift[0], v[1] - shift[1]) for v in Q)))
    return sorted(out)


def expand_sym(D, N):
    """from symmetry-reduced direct enumeration (polygon, weight): all images of each piece under the 8 symmetries
    (deduplicated mod N); returns the list of distinct polygons"""
    seen = {}
    for P, w in D:
        for g in SYMS:
            Q = clean([g(v) for v in P])
            c = centroid_pt(Q)
            shift = (N * floor(c[0] / N), N * floor(c[1] / N))
            key = tuple(sorted((v[0] - shift[0], v[1] - shift[1]) for v in Q))
            seen[key] = Q
    return list(seen.values())


for k in KS:
    N = 5 ** k
    for r in R_LIST:
        t0 = time.time()
        L = lift_components(k, r)
        Ds = direct_components(k, r, sym=True)
        kl = keyset(L, N)
        ks = keyset(expand_sym(Ds, N), N)
        agree = (kl == ks)
        extra = ""
        if k <= 2:
            Df = direct_components(k, r, sym=False)
            extra = f", full direct agrees: {keyset([P for P, w in Df], N) == kl}"
        tot = sum(area(P) for P in L)
        print(f"\n### N = {N}, r = {r} (s = {F(1,2) - r}): {len(L)} components, total area {tot} = {float(tot):.8f}; "
              f"lift = direct(sym) componentwise: {agree}{extra}")
        orbits = defaultdict(list)
        for P in L:
            orbits[orbit_key(P, N)].append(P)
        rows = []
        for key, Ps in orbits.items():
            P = clean(Ps[0])
            kap, opt, how = component_kappa_fast(P, k)
            name, TP, dmax, dmin = classify(P, k)
            ar = area(P)
            assert all(area(clean(Q)) == ar for Q in Ps)
            c = centroid_pt(P)
            assert kappa_point(c, k) >= r
            rows.append((kap, len(Ps), c, ar, name, TP, dmax, opt, len(P)))
        rows.sort(key=lambda t: (-t[0], t[4]))
        print("| kappa | orbit | point c (exact) | area | nearest type pt | max dist | #vert | kappa attained at |")
        print("|---|---|---|---|---|---|---|---|")
        for kap, n, c, ar, name, TP, dmax, opt, nv in rows:
            tname = "c*" if name == "c" else "q"
            print(f"| {kap} = {float(kap):.6f} | {n} | ({c[0]}, {c[1]}) ~ ({float(c[0]):.4f}, {float(c[1]):.4f}) | "
                  f"{ar} ~ {float(ar):.3e} | {tname} ({TP[0]}, {TP[1]}) | {math.sqrt(dmax):.5f} | {nv} | "
                  f"({opt[0]}, {opt[1]}) |")
        print(f"(time {time.time() - t0:.1f}s)")
        sys.stdout.flush()
