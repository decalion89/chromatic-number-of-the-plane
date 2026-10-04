"""Dual certificates for the eight index vectors of family 7 in the window G(2,1) (modulus 325, r0 = 249/1000):
for each, three of its 60 margins and rationals lam >= 0, sum lam = 1, whose combination is the constant 2/7 (the
linear parts cancel).  So no point with these strip indices has all margins > 2/7; with the 7-point (all margins
>= 2/7) the largest least margin of each of these components is exactly 2/7.  Writes cert_K2M1_seven.txt.

Method: at the 7-point c* listed in cert_K2M1.txt, take the margins equal to 2/7 and search the triples (and pairs)
whose gradients have 0 in their convex hull, exactly (fractions)."""
from fractions import Fraction as Fr
from itertools import combinations

def gdiv(a, b):
    n = b[0] ** 2 + b[1] ** 2
    return (Fr(a[0] * b[0] + a[1] * b[1], n), Fr(a[1] * b[0] - a[0] * b[1], n))
def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def gpow(z, e):
    w = (Fr(1), Fr(0))
    if e < 0:
        z, e = gdiv((1, 0), z), -e
    for _ in range(e):
        w = gmul(w, z)
    return w
RHO, SIG = gdiv((2, 1), (2, -1)), gdiv((3, 2), (3, -2))
lines = open("cert_K2M1.txt").read().splitlines()
keys = [tuple(x.strip("()").split(",")) for x in lines[1].split(":", 1)[1].split()]
keys = [(int(j), int(l), e) for (j, l, e) in keys]
def fun(j, l, e):                     # f(x + iy) = a x + b y  for Re or Im of conj(c) rho^j sigma^l
    g = gmul(gpow(RHO, j), gpow(SIG, l))
    return (g[0], g[1]) if e == "Re" else (g[1], -g[0])
F = [fun(*k) for k in keys]
out = [lines[0].replace("29 components", "the 8 components of family 7") + "; dual certificates with value 2/7", lines[1]]
for ln in lines[2:]:
    if "SEVEN" not in ln:
        continue
    idx = [int(t) for t in ln.split("|")[0].split("indices")[1].split()]
    x, y = (Fr(t) for t in ln.split("SEVEN")[1].split())
    marg = []                         # (key, sign, gradient, constant): margin = gradient . c + constant
    for k, (a, b), n in zip(keys, F, idx):
        marg.append((k, "+", (a, b), Fr(-n)))
        marg.append((k, "-", (-a, -b), Fr(n + 1)))
    val = [g[0] * x + g[1] * y + c0 for (_, _, g, c0) in marg]
    assert min(val) == Fr(2, 7), min(val)
    act = [m for m, v in zip(marg, val) if v == Fr(2, 7)]
    cert = None
    for r in (2, 3):
        for tri in combinations(act, r):
            G = [m[2] for m in tri]
            if r == 2:
                # lam g1 + (1 - lam) g2 = 0
                (a1, b1), (a2, b2) = G
                if a1 * b2 - a2 * b1 != 0: continue
                d = (a1 - a2) if a1 != a2 else (b1 - b2)
                if d == 0: continue
                lam = (-a2 / d) if a1 != a2 else (-b2 / d)
                lams = [lam, 1 - lam]
                if not all(l >= 0 for l in lams) or lam * a1 + (1 - lam) * a2 != 0 or lam * b1 + (1 - lam) * b2 != 0: continue
            else:
                (a1, b1), (a2, b2), (a3, b3) = G
                # solve l1 g1 + l2 g2 + l3 g3 = 0, l1 + l2 + l3 = 1
                det = a1 * (b2 - b3) - a2 * (b1 - b3) + a3 * (b1 - b2)
                if det == 0: continue
                l1 = (a2 * b3 - a3 * b2) / det; l2 = (a3 * b1 - a1 * b3) / det; l3 = (a1 * b2 - a2 * b1) / det
                lams = [l1, l2, l3]
                if not all(l >= 0 for l in lams): continue
            const = sum(l * m[3] for l, m in zip(lams, tri))
            gx = sum(l * m[2][0] for l, m in zip(lams, tri)); gy = sum(l * m[2][1] for l, m in zip(lams, tri))
            if gx == 0 and gy == 0 and const == Fr(2, 7) and sum(lams) == 1:
                cert = list(zip(tri, lams)); break
        if cert: break
    assert cert, "no certificate found"
    c = " ; ".join(f"{m[0][0]} {m[0][1]} {m[0][2]} {m[1]} {l}" for m, l in cert)
    out.append(ln.split("|")[0] + f"| SEVEN {x} {y} | kappa = 2/7 ; cert {c}")
open("cert_K2M1_seven.txt", "w").write("\n".join(out) + "\n")
print(f"wrote cert_K2M1_seven.txt: {len(out) - 2} certificates")
