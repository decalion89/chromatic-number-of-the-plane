"""overlap.py d STATE.json [K]: rotations u = w conj(v) for w, v in U_D, ranked by how many points of the sample (the
2000 points of G nearest the origin) they map into G; prints the best K with t = Im(u)/(1 - Re(u)), so that
u = (t + i)/(t - i) (rotunion.py takes t)."""
import sys, json
from fractions import Fraction as Fr
import numpy as np
d = int(sys.argv[1]); st = json.load(open(sys.argv[2])); K = int(sys.argv[3]) if len(sys.argv) > 3 else 10
D = st["D"]; U = np.array(st["U"], dtype=np.int64); P = np.array(st["P"], dtype=np.int64)
Pset = {tuple(p) for p in P.tolist()}
r = (P[:, 0] + P[:, 1] * d ** 0.5) ** 2 + (P[:, 2] + P[:, 3] * d ** 0.5) ** 2
S = P[np.argsort(r)[:2000]]
def cmul(x, y):
    """(a + b s) + i (c + e s) times the same, s = sqrt d, integer coefficients"""
    a1, b1, c1, e1 = x; a2, b2, c2, e2 = y
    re0 = a1 * a2 + d * b1 * b2 - (c1 * c2 + d * e1 * e2); re1 = a1 * b2 + b1 * a2 - (c1 * e2 + e1 * c2)
    im0 = a1 * c2 + d * b1 * e2 + c1 * a2 + d * e1 * b2; im1 = a1 * e2 + b1 * c2 + c1 * b2 + e1 * a2
    return (re0, re1, im0, im1)
seen, res = set(), []
Ul = [tuple(u) for u in U.tolist()]
for w in Ul:
    for v in Ul:
        vb = (v[0], v[1], -v[2], -v[3])
        u = cmul(w, vb)                                  # over D^2
        if u in seen: continue
        seen.add(u)
        if u == (D * D, 0, 0, 0): continue
        # u * p over D^3 for the sample; in G iff divisible by D^2
        a1, b1, c1, e1 = u
        a2, b2, c2, e2 = S[:, 0], S[:, 1], S[:, 2], S[:, 3]
        re0 = a1 * a2 + d * b1 * b2 - (c1 * c2 + d * e1 * e2); re1 = a1 * b2 + b1 * a2 - (c1 * e2 + e1 * c2)
        im0 = a1 * c2 + d * b1 * e2 + c1 * a2 + d * e1 * b2; im1 = a1 * e2 + b1 * c2 + c1 * b2 + e1 * a2
        Q = np.stack([re0, re1, im0, im1], axis=1)
        ok = np.all(Q % (D * D) == 0, axis=1)
        hits = sum(1 for q in (Q[ok] // (D * D)).tolist() if tuple(q) in Pset)
        res.append((hits, u))
res.sort(reverse=True)
print(f"{len(res)} distinct rotations; best overlaps on the 2000-point sample:")
for hits, u in res[:K]:
    x = (Fr(u[0], D * D), Fr(u[1], D * D)); y = (Fr(u[2], D * D), Fr(u[3], D * D))
    # t = y / (1 - x) in K
    a, b = 1 - x[0], -x[1]; n = a * a - d * b * b
    inv = (a / n, -b / n)
    t = (y[0] * inv[0] + d * y[1] * inv[1], y[0] * inv[1] + y[1] * inv[0])
    print(f"  overlap {hits}: u = {u}/{D}^2   t = {t[0]} {t[1]}")
