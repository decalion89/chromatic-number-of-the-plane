"""Independent exact computation (own code, cell-by-cell polygon clipping with Fractions) of the components of
S^theta(2,1) = {c : Re(conj(c) gamma) in [theta, 1-theta] + Z for all gamma in G(2,1)}, modulo 325 Z[i].

Usage: python3 window21_components.py THETA [a_start a_end]   (THETA like 2/7; cells a in [a_start, a_end))
Writes one line per surviving component: cell, index vector (30 strip indices), vertices, contained type points.

The 30 functionals: (j,l) in [-2..2]x[-1..1], Re and Im of conj(c) rho^j sigma^l, where for c = x+iy and
g = 325 rho^j sigma^l = g0 + i g1 in Z[i]:  Re(conj(c) g)/325 = (g0 x + g1 y)/325,  Im(conj(c) g)/325 = (g1 x - g0 y)/325.
The conditions at -gamma and +-i gamma are equivalent to these (the interval is symmetric), so these 30 suffice.
(0,0,Re) is x and (0,0,Im) is -y, so the cells are [a+theta, a+1-theta] x [b+theta, b+1-theta], 0 <= a, b < 325.
A connected component lies in one cell (x and y stay in one strip), and is the convex polygon of one index vector.
"""
import sys
from fractions import Fraction as Fr

def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def gpow(z, e):
    r = (1, 0)
    for _ in range(e):
        r = gmul(r, z)
    return r

def g325(j, l):
    # 325 rho^j sigma^l with rho = (2+i)/(2-i) = (2+i)^2/5, sigma = (3+2i)/(3-2i) = (3+2i)^2/13
    # 325 rho^j sigma^l = (2+i)^(2+j) (2-i)^(2-j) (3+2i)^(1+l) (3-2i)^(1-l)
    g = (1, 0)
    for z, e in (((2, 1), 2 + j), ((2, -1), 2 - j), ((3, 2), 1 + l), ((3, -2), 1 - l)):
        g = gmul(g, gpow(z, e))
    return g

# sanity: these are 325 times the rotations (3+4i)/5 and (5+12i)/13
assert g325(1, 0) == (195, 260) and g325(0, 1) == (125, 300) and g325(0, 0) == (325, 0)

KEYS = [(j, l, p) for j in range(-2, 3) for l in range(-1, 2) for p in ("Re", "Im")]
FORMS = {}
for (j, l, p) in KEYS:
    g0, g1 = g325(j, l)
    FORMS[(j, l, p)] = (g0, g1) if p == "Re" else (g1, -g0)
assert FORMS[(0, 0, "Re")] == (325, 0) and FORMS[(0, 0, "Im")] == (0, -325)
ORDER = [k for k in KEYS if k[:2] != (0, 0)]
# order: rotations far from 1 first (they cut the axis-parallel cells most)
ORDER.sort(key=lambda k: (-abs(k[0]) - abs(k[1]), k))

def clip(poly, A, B, C):
    """intersect convex polygon (list of Fraction points, possibly degenerate) with A x + B y <= C"""
    out = []
    n = len(poly)
    vals = [A * P[0] + B * P[1] - C for P in poly]
    for i in range(n):
        P, vP = poly[i], vals[i]
        Q, vQ = poly[(i + 1) % n], vals[(i + 1) % n]
        if vP <= 0:
            out.append(P)
        if (vP < 0 < vQ) or (vQ < 0 < vP):
            t = vP / (vP - vQ)
            out.append((P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1])))
    # dedupe
    res = []
    for P in out:
        if not res or res[-1] != P:
            res.append(P)
    while len(res) > 1 and res[0] == res[-1]:
        res.pop()
    return res

def ceil_fr(q):
    return -((-q.numerator) // q.denominator)

def floor_fr(q):
    return q.numerator // q.denominator

def search(theta, a_range, out):
    one_m = 1 - theta
    found = 0
    for a in a_range:
        for b in range(325):
            cell = [(a + theta, b + theta), (a + one_m, b + theta), (a + one_m, b + one_m), (a + theta, b + one_m)]
            idx0 = {(0, 0, "Re"): a, (0, 0, "Im"): -b - 1}
            stack = [(0, cell, idx0)]
            while stack:
                d, poly, idx = stack.pop()
                if d == len(ORDER):
                    found += 1
                    out.append((a, b, idx, poly))
                    continue
                k = ORDER[d]
                g0, g1 = FORMS[k]
                vals = [(g0 * P[0] + g1 * P[1]) / 325 for P in poly]
                fmin, fmax = min(vals), max(vals)
                nlo = ceil_fr(fmin - one_m)
                nhi = floor_fr(fmax - theta)
                for n in range(nlo, nhi + 1):
                    # n + theta <= f <= n + 1 - theta, f = (g0 x + g1 y)/325
                    p1 = clip(poly, -g0, -g1, -325 * (n + theta))
                    if not p1:
                        continue
                    p2 = clip(p1, g0, g1, 325 * (n + one_m))
                    if not p2:
                        continue
                    idx2 = dict(idx)
                    idx2[k] = n
                    stack.append((d + 1, p2, idx2))
    return found

def type_points():
    pts = []
    pts.append(("c", (Fr(325, 2), Fr(325, 2))))
    for al in (1, 2):
        for be in (1, 2):
            pts.append((f"q({al},{be})", (Fr(325 * al, 3), Fr(325 * be, 3))))
    for a in range(7):
        for b in range(7):
            pts.append((f"7?({a},{b})", (Fr(325 * a, 7), Fr(325 * b, 7))))
    return pts

def in_poly(pt, poly, theta):
    # check the point satisfies all constraints of the component: easiest is to test floors & margins directly
    return None

if __name__ == "__main__":
    theta = Fr(sys.argv[1])
    a0, a1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, 325)
    out = []
    import time
    t = time.time()
    found = search(theta, range(a0, a1), out)
    print(f"# theta={theta} cells a in [{a0},{a1}) : {found} components  [{time.time() - t:.1f}s]", flush=True)
    for a, b, idx, poly in out:
        vec = [idx[k] for k in KEYS]
        print("COMP", a, b, " ".join(map(str, vec)), "|", " ".join(f"{P[0]},{P[1]}" for P in poly))
