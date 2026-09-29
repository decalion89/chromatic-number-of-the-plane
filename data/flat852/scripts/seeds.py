import os
"""seeds.py -- seed point sets (integer 12-vectors, scaled by 7, power basis of zeta21)."""
import json, cmath, math
import numpy as np
from flat import *

SC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the scratchpad directory
KS = [k for k in range(1, 21) if math.gcd(k, 21) == 1]          # the 12 embeddings z -> exp(2 pi i k/21)
EMB = np.array([[cmath.exp(2j * math.pi * j * k / 21) for k in KS] for j in range(N)])   # (12,12)


def from_embeddings(vals):
    """integer vector c with sum_j c_j exp(2 pi i j k/21) = 7*vals[k] for the 12 embeddings (rounded)"""
    A = np.concatenate([EMB.T.real, EMB.T.imag])
    b = 7 * np.concatenate([np.real(vals), np.imag(vals)])
    c, *_ = np.linalg.lstsq(A, b, rcond=None)
    ci = np.rint(c).astype(np.int64)
    assert np.max(np.abs(c - ci)) < 1e-6, np.max(np.abs(c - ci))
    return ci


def heptagon():
    """the 21 points of Haugland's H (centre at the origin), exact; checked against the paper's formulas"""
    pts = {}
    for name, fn in (("P", lambda s: -1 / (s(3) - 1 / s(3))),                  # z = zeta42 in hept.py
                     ("Q", lambda s: s(14) / (s(6) - 1 / s(6))),
                     ("R", lambda s: s(28) / (s(12) - 1 / s(12)))):
        # sigma_k(zeta42) : zeta42 = -zeta21^11 -> -exp(2 pi i 11 k/21)
        vals = []
        for k in KS:
            z42 = -cmath.exp(2j * math.pi * 11 * k / 21)
            s = lambda e, z42=z42: z42 ** e
            vals.append(fn(s))
        base = from_embeddings(np.array(vals))
        for j in range(7):
            # multiply by zeta7^j = zeta21^(3j)
            pts[f"{name}{j}"] = np.array(mul(list(base), zpow(3 * j)), dtype=np.int64)
    # check against the paper's numeric formulas
    al, be, ga = 1 / math.sin(2 * math.pi / 7), 1 / math.sin(4 * math.pi / 7), 1 / math.sin(8 * math.pi / 7)
    for j in range(7):
        t = 2 * math.pi * j / 7
        for nm, val in (("P", (-ga / 2) * cmath.exp(1j * (t + math.pi / 2))),
                        ("Q", (al / 2) * cmath.exp(1j * (t + math.pi / 6))),
                        ("R", (be / 2) * cmath.exp(1j * (t + 5 * math.pi / 6)))):
            assert abs(cval(pts[f"{nm}{j}"]) - val) < 1e-9, (nm, j)
    return pts


def g1_points():
    """Haugland's G1 as rebuilt in record_hept (basis zeta42 there), converted to the zeta21 basis"""
    d = json.load(open(f"{SC}/record_hept/G1.json"))
    out = []
    for v in d["vertices"]:
        c = v["exact7"] if "exact7" in v else v["exact"]
        acc = [0] * N
        for k, x in enumerate(c):
            if x:
                t = [(-1) ** k * x * y for y in zpow(11 * k)]
                acc = [a + b for a, b in zip(acc, t)]
        out.append(acc)
        assert abs(cval(acc) - complex(float(v["x"]), float(v["y"]))) < 1e-9
    return np.array(out, dtype=np.int64), d["A"], d["B"]


if __name__ == "__main__":
    D, CD, U = directions()
    H = heptagon()
    P = np.array(list(H.values()))
    E, J = build_edges(P, U)
    Uset = set(D)
    print("H: 21 points, U-edges:", len(E), " all in D:", all(tuple(P[b] - P[a]) in Uset or tuple(P[a] - P[b]) in Uset for a, b in E))
    inO = [k for k, v in H.items() if np.all(v % 7 == 0)]
    print("points of H in Z[zeta21]:", inO)
    C = conj_rows(P)
    Q = unique_rows(np.concatenate([P, C]))
    E2, _ = build_edges(Q, U)
    print("H u conj(H):", len(Q), "points,", len(E2), "U-edges")
    G, A, B = g1_points()
    EG, JG = build_edges(G, U)
    print("G1:", len(G), "points,", len(EG), "U-edges; all in D:",
          all(tuple(G[b] - G[a]) in Uset for a, b in EG))
    Q = unique_rows(np.concatenate([G, conj_rows(G)]))
    E3, _ = build_edges(Q, U)
    print("G1 u conj(G1):", len(Q), "points,", len(E3), "U-edges")
