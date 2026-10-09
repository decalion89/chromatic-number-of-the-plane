"""Exact a-posteriori check that every computed polygon is a whole connected component of S_N^(r): all vertices have
kappa >= r (kappa computed from the definition, min over the 4(2k+1) rotations), and every edge midpoint of a
2-dimensional polygon has kappa exactly r (so each edge lies on a supporting line of the component, hence the polygon
equals the component).  Degenerate components (points, segments) are checked for kappa >= r only."""
from snr import *
import sys
R_LIST = [F(3001, 10000), F(302, 1000), F(305, 1000), F(31, 100), F(32, 100), F(1, 3), F(17, 56), F(3, 10)]
for k in (1, 2, 3):
    for r in R_LIST:
        L = [clean(P) for P in lift_components(k, r)]
        ok = True; nd = 0
        for P in L:
            ok &= all(kappa_point(v, k) >= r for v in P)
            if len(P) >= 3:
                for i in range(len(P)):
                    a, b = P[i], P[(i + 1) % len(P)]
                    m = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
                    ok &= kappa_point(m, k) == r
            else:
                nd += 1
        print(f"N={5**k} r={r}: {len(L)} components ({nd} degenerate): vertices in S and edges on the boundary: {ok}")
        sys.stdout.flush()
