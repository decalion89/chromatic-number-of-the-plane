"""Does a colouring factor through a place?  Residues of points at unramified places vs colours.

For each prime p (unramified in the field) and each place (choice of square roots in Q_{p^2}),
every point z gets its residue g(z) in F_{p^2} relative to its coset of the p-adic integers
(units digits, as in hn/adelic.py).  A colouring that factors through that place has each
residue class monochromatic: purity = fraction of points whose colour is the majority colour of
their class.  Purity near 1 flags a hidden periodic colouring; near 1/5 means none.

usage: python3 placediag.py <checkpoint.json with points and colouring> <p1,p2,...>"""
import sys, os, json, itertools
from collections import Counter, defaultdict
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from hn.field import Field
from hn.geometry import Point
from residuegate import sqrt_in_Qp2, mul

N = 12


def units_digit(num, den, p, m):
    # r = num/den; (r - frac_p(r)) mod p
    k = 0
    while den % p == 0:
        den //= p; k += 1
    mm = p ** (k + 1)
    return ((num * pow(den, -1, mm)) % mm) // (p ** k)


def point_residues(P, F, p, r, signs):
    m = p ** N
    roots = {}
    for g, sg in zip(list(F.gens) + [-1], signs):
        a, b = sqrt_in_Qp2(g % m if g > 0 else (m + g), p, r, N)
        roots[g] = ((sg * a) % m, (sg * b) % m)
    B = []
    for i in range(F.dim):
        v = (1, 0)
        for j, g in enumerate(F.gens):
            if i >> j & 1:
                v = mul(v, roots[g], r, m)
        B.append(v)
    I = roots[-1]
    out = []
    for z in P:
        # z = sum c_i B_i + i * sum d_i B_i, with rational c, d: accumulate A + B sqrt r with a
        # common denominator, then take units digits
        den = 1
        terms = []
        for part, extra in ((z.x.c, (1, 0)), (z.y.c, I)):
            for i, c in enumerate(part):
                c = Fr(c)
                if c:
                    t = mul(B[i], extra, r, m)
                    terms.append((c, t))
                    den = den * c.denominator // __import__("math").gcd(den, c.denominator)
        A = sum(int(c * den) * t[0] for c, t in terms)
        Bv = sum(int(c * den) * t[1] for c, t in terms)
        out.append((units_digit(A, den, p, m), units_digit(Bv, den, p, m)))
    return out


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    F = Field(tuple(d["field_generators"]))
    mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
    col = d["colouring"]
    P = [mk(xy) for xy in d["points"]][:len(col)]
    print(f"{os.path.basename(sys.argv[1])}: {len(P)} coloured points", flush=True)
    for p in [int(x) for x in sys.argv[2].split(",")]:
        if any(g % p == 0 for g in F.gens) or p == 2:
            print(f"  p = {p}: ramified or 2, skipped"); continue
        r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        best = 0.0
        for signs in itertools.product((1, -1), repeat=len(F.gens) + 1):
            res = point_residues(P, F, p, r, signs)
            cls = defaultdict(Counter)
            for g, c in zip(res, col):
                cls[g][c] += 1
            purity = sum(max(cn.values()) for cn in cls.values()) / len(P)
            best = max(best, purity)
        print(f"  p = {p}: best purity over the places above p: {best:.3f}  (classes {len(cls)})", flush=True)
