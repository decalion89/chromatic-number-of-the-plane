"""Consistency check of Theorem A: for r slightly above 17/56 the lift enumeration finds exactly the 5 main components
at every level k <= K; also checks the shapes: the c*-component equals c* + P_k^(s) (exact vertex sets) and the
q-components have radius <= sqrt2 * (1/3 - r)."""
from snr import *
from radii import polygon_from_halfplanes
from lemma_lp import lin
import sys
r = F(sys.argv[1]); K = int(sys.argv[2])
s = F(1, 2) - r; eps = F(1, 3) - r
lo, hi = r, 1 - r
comps = [[(lo, lo), (hi, lo), (hi, hi), (lo, hi)]]
for kk in range(1, K + 1):
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
    cs, qs = type_points(kk)
    # P_k^(s)
    hp = []
    for j in range(-kk, kk + 1):
        for (a, b) in lin(j):
            hp.append((a, b, -s)); hp.append((-a, -b, -s))
    Pk = polygon_from_halfplanes(hp)
    ok_c = ok_q = 0
    for P in comps:
        P = clean(P)
        name, TP, dmax, dmin = classify(P, kk)
        if name == "c":
            shifted = sorted(((v[0] - TP[0]) - N * floor((v[0] - TP[0]) / N + F(1, 2)), (v[1] - TP[1]) - N * floor((v[1] - TP[1]) / N + F(1, 2))) for v in P)
            ok_c += shifted == sorted(Pk)
        else:
            ok_q += dmax <= 2 * eps * eps
    print(f"level {kk}: {len(comps)} components; c*-component == c* + P_k^(s): {ok_c == 1}; q-components with radius <= sqrt2*eps: {ok_q}")
