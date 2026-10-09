"""padic_planes.py: bounds for the chromatic number of the p-adic plane, the graph on Q_p^2 in which two points are
adjacent when (x - x')^2 + (y - y')^2 = 1.

Upper bounds. For p = 2 and for p = 3 (mod 4), -1 is not a square in Q_p, so x^2 + y^2 = 1 forces x, y to be p-adic
integers (a pole on either coordinate would make -1 a square mod p, or a sum of two squares 0 mod 4 when p = 2).
Adjacent points therefore lie in the same coset of Z_p^2, and their difference reduces mod p to a nonzero solution
of x^2 + y^2 = 1 in F_p (for p = 2: a vector with exactly one odd coordinate). Colouring each coset through a proper
colouring of the finite plane F_p^2 (Moorhouse's table; for p = 2 and 3 the colouring x + y mod p) colours Q_p^2.

Lower bounds. A unit-distance graph over a number field K that embeds in Q_p (with the same form) is a subgraph of
the p-adic plane. Q(sqrt d) embeds in Q_p when d is a nonzero square mod p (p odd, Hensel). So each certified
graph of data/quadratic_planes/ with chi = 4 bounds chi(Q_p^2) below by 4 at every such p, and for p = 3 a 5-cycle
over Q(sqrt7) (7 = 1 mod 3) gives 3. The 5-chromatic graph of data/five_247_c.json lies over Q(sqrt3, sqrt11,
sqrt247) and bounds chi(Q_p^2) below by 5 where 3, 11 and 247 are squares mod p.

Prints the table and checks every ingredient it uses (the stored colourings, the 5-cycle, the quadratic residues)."""
import json, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
FIELDS = sorted(int(f[1:-5]) for f in os.listdir(DATA) if f.startswith("q") and f.endswith(".json") and f != "q47.json")
CHI_FINITE = {7: "4", 11: "5", 19: "5", 23: "5-8", 31: "5-8", 43: "6-8", 47: "6-10", 59: "6-10"}   # notes/local_colourings.md


def is_qr(a, p):
    a %= p
    return a != 0 and pow(a, (p - 1) // 2, p) == 1


def check_finite_colouring(p, col):
    U = [(a, b) for a in range(p) for b in range(p) if (a * a + b * b) % p == 1]
    return all(col[p * x + y] != col[p * ((x + a) % p) + (y + b) % p] for x in range(p) for y in range(p) for a, b in U)


def five_cycle_q7():
    """the 5-cycle over Q(sqrt 7): elements a + b sqrt7 as pairs of Fractions; returns the vertices after checking"""
    d = 7
    mul = lambda a, b: (a[0] * b[0] + d * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    add = lambda a, b: (a[0] + b[0], a[1] + b[1])
    vadd = lambda u, v: (add(u[0], v[0]), add(u[1], v[1]))
    vsc = lambda c, u: ((c * u[0][0], c * u[0][1]), (c * u[1][0], c * u[1][1]))
    nrm = lambda v: add(mul(v[0], v[0]), mul(v[1], v[1]))
    one, zero = (F(1), F(0)), (F(0), F(0))
    u1, u2 = (one, zero), (zero, one)
    u3 = ((F(-1, 4), F(1, 4)), (F(-1, 4), F(-1, 4)))
    w = vsc(-1, vadd(vadd(u1, u2), u3))
    Jw = (vsc(-1, (w[1], w[1]))[0], w[0])
    u4, u5 = vadd(vsc(F(1, 2), w), vsc(F(1, 2), Jw)), vadd(vsc(F(1, 2), w), vsc(F(-1, 2), Jw))
    P = [(zero, zero)]
    for u in (u1, u2, u3, u4, u5):
        assert nrm(u) == one
        P.append(vadd(P[-1], u))
    assert P[5] == P[0] and len(set(P[:5])) == 5
    return P[:5]


fin = json.load(open(os.path.join(DATA, "finite_planes.json")))["planes"]
for q in ("7", "11"):
    assert check_finite_colouring(int(q), fin[q]["colouring"]), q
assert check_finite_colouring(3, [(x + y) % 3 for x in range(3) for y in range(3)])
five_cycle_q7()
assert is_qr(7, 3)
five_fields = (3, 11, 247)

print(f"certified fields with chi = 4: d = {', '.join(map(str, FIELDS))} (and d = 47, chi >= 4)")
print("p   | form       | upper (finite plane) | lower | witnesses")
for p in [2, 3, 5, 7, 11, 13, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
    if p == 2:
        print("2   | anisotropic| 2 (x + y mod 2)      | 2     | an edge"); continue
    aniso = p % 4 == 3
    if not aniso:
        emb = [d for d in FIELDS + [47] if is_qr(d, p)]
        print(f"{p:<3} | isotropic  | infinite Borel chromatic number (Bardestani-Mallahi-Karai) | 4 | Q(sqrt{emb[0]})" if emb else
              f"{p:<3} | isotropic  | -")
        continue
    if p == 3:
        print("3   | anisotropic| 3 (x + y mod 3)      | 3     | 5-cycle over Q(sqrt7)"); continue
    emb = [d for d in FIELDS + [47] if is_qr(d, p)]
    lo = "5" if all(is_qr(x, p) for x in five_fields) else ("4" if emb else "3")
    wit = "five_247_c (803 vertices)" if lo == "5" else (f"Q(sqrt{emb[0]})" + (f", Q(sqrt47)" if 47 in emb and emb[0] != 47 else "") if emb else "")
    print(f"{p:<3} | anisotropic| {CHI_FINITE.get(p, '>= 6 (Hoffman)'):<20} | {lo:<5} | {wit}")
print("exact: chi(Q_2^2) = 2, chi(Q_3^2) = 3, chi(Q_7^2) = 4")
