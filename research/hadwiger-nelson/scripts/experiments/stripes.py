"""Stripe-twisted coset colourings:  c(p) = psi(p) + k * floor(<p, v> / w)  (mod 5).

With psi a coset colouring (psi(u) != 0 on every unit vector) and w >= 1, a
unit step changes the stripe index by 0 or +-1, so c is proper exactly when
every unit vector u with psi(u) = k has <u, v> >= 0 -- the class psi = k fits
in a closed half-plane.  Such a colouring is NOT constant on cosets: across a
stripe boundary it shifts by k, so it splits pairs u, u + 5e (and can keep
2e pairs alike).  One colouring of the whole module, refuting forcing on
every graph drawn in it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, itertools, math, time
from fractions import Fraction as Fr
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(HN_DIR + "/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = HN_DIR
NAME = sys.argv[1]
t0 = time.time()
if NAME == "ei":
    sols = [(a,b,c,d) for a in range(-7,8) for b in range(-4,5) for c in range(-12,13) for d in range(-3,4)
            if 3*a*a + 11*b*b + c*c + 33*d*d == 144 and a*b == -c*d]
    E = list({max(v, tuple(-x for x in v)) for v in sols})
    ang = {v: math.atan2((v[2] + v[3] * math.sqrt(33)) / 12, (v[0] * math.sqrt(3) + v[1] * math.sqrt(11)) / 12) for v in E}
else:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P)
    # angles of the edge directions, keyed by the same integer vectors edge_vectors produces
    raw = {}
    for a, b in g.edges():
        dx = g.vertices[b].x - g.vertices[a].x; dy = g.vertices[b].y - g.vertices[a].y
        raw[tuple(dx.c) + tuple(dy.c)] = math.atan2(float(dy), float(dx))
    E = edge_vectors(g)
    den = 1
    from math import gcd
    for v in raw:
        for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
    ang = {}
    for v, th in raw.items():
        iv = tuple(int(Fr(q) * den) for q in v)
        key = max(iv, tuple(-x for x in iv))
        ang[key] = th if key == iv else math.atan2(-math.sin(th), -math.cos(th))
B = echelon(E); C = np.array([coords(B, v) for v in E]); r = len(B)
TH = np.array([ang[v] for v in E])
units = [(C[k], TH[k]) for k in range(len(E))] + [(-C[k], math.atan2(-math.sin(TH[k]), -math.cos(TH[k]))) for k in range(len(E))]
PS = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
ok = np.ones(len(PS), bool)
for c in C: ok &= (PS @ c) % 5 != 0
ADM = PS[ok]
print(f"{NAME}: {len(E)} directions, rank {r}, {len(ADM)} admissible coset colourings", flush=True)
UC = np.array([u for u, th in units]); UT = np.array([th for u, th in units])
def fits_half_plane(angles):
    if len(angles) == 0: return True, 0.0
    a = np.sort(np.mod(angles, 2 * math.pi))
    gaps = np.diff(np.concatenate([a, [a[0] + 2 * math.pi]]))
    i = int(np.argmax(gaps))
    if gaps[i] >= math.pi - 1e-12:           # all within a closed semicircle
        mid = (a[i] + gaps[i] / 2) % (2 * math.pi)   # centre of the empty arc
        return True, (mid + math.pi) % (2 * math.pi)  # v points away from the gap
    return False, None
found = []
vals = (ADM @ UC.T) % 5
for i, psi in enumerate(ADM):
    for k in range(1, 5):
        cls = UT[vals[i] == k]
        okk, v = fits_half_plane(cls)
        if okk:
            found.append((tuple(int(x) for x in psi), k, v, len(cls)))
            break
print(f"  stripe-twisted colourings: {len(found)} of {len(ADM)} admissible psi admit one", flush=True)
if found:
    psi, k, v, m = found[0]
    print(f"  e.g. psi={psi}, k={k}, class size {m}, stripe normal at {math.degrees(v):.2f} degrees", flush=True)
