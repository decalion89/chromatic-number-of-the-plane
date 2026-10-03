"""structure_lemma.py [KMAX] [KPOLY]: exact checks of every computational fact used in the proof that
chi(Q(sqrt d)^2) >= 4 for d = 11 (mod 12) (notes/four_colours_11_mod_12.md, Proposition 1).

Notation (all exact, Gaussian rationals as pairs of Fractions):
  N = 5^k, rho = (3+4i)/5, G_N = rational unit vectors with denominator dividing N,
  S_N = {c in C : Re(conj(c) g) in [1/3, 2/3] + Z for every g in G_N},
  c* = N(1+i)/2,  Q_N = {(N/3)(a+bi) : a, b in {1, 2}},  Sq0 = [-1/6, 1/6]^2,
  P_k = {x : conj(x) rho^j in Sq0 for |j| <= k},  Lambda_+ = ((2+i)/5) Z[i].
Checks:
  (1) G_N = {e rho^j : e in {1, i, -1, -i}, |j| <= k}, |G_N| = 4(2k+1), for k <= KMAX;
  (2) Re and Im of conj(c*) g lie in 1/2 + Z for every g in G_N (c* is the bipartite character);
  (3) Q_N is contained in S_N;
  (4) the vertices of P_1 and max |x|^2 over P_1 = 10/324 = (sqrt10/18)^2 < 18/324 = (sqrt2/6)^2;
  (5) the four nonzero classes of Lambda_+ / Z[i] are at squared distance exactly 1/18 from Sq0 + Z[i]
      (and no point of Lambda_+ outside Z[i] lies in [-1/3, 1/3]^2 + Z[i]); the same for Lambda_- = conj;
  (6) for k <= KPOLY: S_{5^k} computed by exact polygon clipping is the polygon c* + P_k plus the 4 points Q_N
      (independent of the induction in the note)."""
import sys
from fractions import Fraction as F
from math import isqrt, gcd

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
KPOLY = int(sys.argv[2]) if len(sys.argv) > 2 else 2
LO, HI = F(1, 3), F(2, 3)


def mul(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def conj(a): return (a[0], -a[1])
def frac(x): return x - (x.numerator // x.denominator)


RHO = (F(3, 5), F(4, 5))
UNITS = [(F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1))]


def rho_pow(j):
    r = RHO if j >= 0 else conj(RHO)
    z = (F(1), F(0))
    for _ in range(abs(j)): z = mul(z, r)
    return z


def rational_units(N):
    out = set()
    for c in range(1, N + 1):
        if N % c: continue
        for a in range(-c, c + 1):
            b2 = c * c - a * a; b = isqrt(b2)
            if b * b == b2:
                for bb in {b, -b}:
                    if gcd(gcd(abs(a), abs(bb)), c) == 1:
                        out.add((F(a, c), F(bb, c)))
    return out


ok = True
def check(cond, msg):
    global ok
    print(("ok   " if cond else "FAIL ") + msg)
    ok = ok and cond


for k in range(0, KMAX + 1):
    N = 5 ** k
    G = rational_units(N)
    want = {mul(e, rho_pow(j)) for e in UNITS for j in range(-k, k + 1)}
    check(G == want and len(G) == 4 * (2 * k + 1), f"(1) k={k}: G_N = e rho^j, |G_N| = {len(G)}")
    cs = (F(N, 2), F(N, 2))
    check(all(frac(z[0]) == F(1, 2) and frac(z[1]) == F(1, 2) for z in (mul(conj(cs), g) for g in G)),
          f"(2) k={k}: conj(c*) g in (1+i)/2 + Z[i] for all g in G_N")
    good = True
    for a in (1, 2):
        for b in (1, 2):
            q = (F(N * a, 3), F(N * b, 3))
            for g in G:
                z = mul(conj(q), g)
                if not (LO <= frac(z[0]) <= HI and LO <= frac(z[1]) <= HI): good = False
    check(good, f"(3) k={k}: Q_N is contained in S_N")

# (4) P_1 = {x : |<x, n>| <= 1/6 for the six normals n}, conj(x) rho^j in Sq0 <=> x in rho^j Sq0
normals = []
for j in (-1, 0, 1):
    r = rho_pow(j)
    normals += [r, mul((F(0), F(1)), r)]
H = F(1, 6)
lines = [(n, s) for n in normals for s in (1, -1)]            # <x, n> = s H
verts = set()
for i in range(len(lines)):
    for j in range(i + 1, len(lines)):
        (n1, s1), (n2, s2) = lines[i], lines[j]
        det = n1[0] * n2[1] - n1[1] * n2[0]
        if det == 0: continue
        x = ((s1 * H) * n2[1] - (s2 * H) * n1[1]) / det
        y = (n1[0] * (s2 * H) - n2[0] * (s1 * H)) / det
        if all(abs(x * n[0] + y * n[1]) <= H for n in normals): verts.add((x, y))
r2 = max(v[0] ** 2 + v[1] ** 2 for v in verts)
check(len(verts) == 12 and r2 == F(10, 324), f"(4) P_1 is a 12-gon, max |x|^2 = {r2} = 10/324 (vertices e.g. {sorted(verts)[0]})")
check(r2 < F(18, 324), "(4) sqrt10/18 < sqrt2/6")


def dist2_to_sq0_lattice(p):
    """squared distance from p to Sq0 + Z[i]"""
    best = None
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            x, y = p[0] - dx, p[1] - dy
            ex = max(F(0), abs(x) - H); ey = max(F(0), abs(y) - H)
            v = ex * ex + ey * ey
            best = v if best is None or v < best else best
    return best


for name, g in (("Lambda_+", (F(2, 5), F(1, 5))), ("Lambda_-", (F(2, 5), F(-1, 5)))):
    classes = [(frac(m * g[0]), frac(m * g[1])) for m in range(1, 5)]
    ds = [dist2_to_sq0_lattice(c) for c in classes]
    check(all(v == F(1, 18) for v in ds), f"(5) {name}: nonzero classes {classes} at squared distance {set(ds)} = 1/18")
    inner = [c for c in classes if (min(frac(c[0]), 1 - frac(c[0])) <= F(1, 3) and min(frac(c[1]), 1 - frac(c[1])) <= F(1, 3))]
    check(not inner, f"(5) {name}: no nonzero class meets [-1/3,1/3]^2 + Z[i]")


# (6) exact polygon clipping of S_{5^k}
def clip(poly, a, b, c, sense):
    out = []; n = len(poly)
    for i in range(n):
        P, Q = poly[i], poly[(i + 1) % n]
        fp = sense * (a * P[0] + b * P[1] - c); fq = sense * (a * Q[0] + b * Q[1] - c)
        if fp >= 0: out.append(P)
        if (fp > 0 and fq < 0) or (fp < 0 and fq > 0):
            t = fp / (fp - fq); out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    res = []
    for p in out:
        if not res or res[-1] != p: res.append(p)
    if len(res) > 1 and res[0] == res[-1]: res.pop()
    return res


def strip(polys, a, b):
    out = []
    for poly in polys:
        vals = [a * p[0] + b * p[1] for p in poly]
        lo, hi = min(vals), max(vals)
        for n in range((lo - HI).__floor__(), (hi - LO).__floor__() + 2):
            q = clip(poly, a, b, n + LO, +1)
            if q: q = clip(q, a, b, n + HI, -1)
            if q: out.append(q)
    return out


def area(poly):
    if len(poly) < 3: return F(0)
    s = F(0)
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2


def in_Pk(x, k):
    for j in range(-k, k + 1):
        z = mul(conj(x), rho_pow(j))
        if abs(z[0]) > H or abs(z[1]) > H: return False
    return True


for k in range(1, KPOLY + 1):
    N = 5 ** k
    polys = [[(n + LO, m + LO), (n + HI, m + LO), (n + HI, m + HI), (n + LO, m + HI)] for n in range(N) for m in range(N)]
    for j in list(range(1, k + 1)) + list(range(-1, -k - 1, -1)):
        A, B = rho_pow(j)
        polys = strip(polys, A, B); polys = strip(polys, B, -A)
    big = [p for p in polys if area(p) > 0]
    pts = [list(dict.fromkeys(p)) for p in polys if area(p) == 0]
    cs = (F(N, 2), F(N, 2))
    Qn = {(F(N * a, 3), F(N * b, 3)) for a in (1, 2) for b in (1, 2)}
    good = len(big) == 1 and all(in_Pk((v[0] - cs[0], v[1] - cs[1]), k) for v in big[0])
    good = good and all(len(p) == 1 for p in pts) and {p[0] for p in pts} == Qn
    # every vertex of P_k, translated by c*, is a vertex of the computed polygon: compare areas with P_k computed directly
    check(good, f"(6) k={k}: S_N = one polygon around c* inside c* + P_k (area {float(area(big[0])) if big else 0:.6f}, "
                f"{len(dict.fromkeys(big[0])) if big else 0} vertices) + the 4 points Q_N")

print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
