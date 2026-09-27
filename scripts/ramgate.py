"""Residue gate at a prime p that may ramify: embed Q(sqrt g_1, ..., sqrt g_n, i) in
L = Q_p(sqrt r, sqrt p) (r a non-residue unit), for every choice of signs, and reduce the unit
vectors mod sqrt p.  Residues lie in F_{p^2} = F_p(sqrt r).  If every unit is integral at the
place, any colouring of Cay(F_{p^2}, residues) colours every graph along these units.

usage: python3 ramgate.py <module.json> <p> [k=5] [TL seconds per SAT]"""
import sys, os, json, itertools
from fractions import Fraction as Fr
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from hn.field import Field
from hn.geometry import Point
from pysat.solvers import Solver

N = 14


def hensel_sqrt(a, p, m):
    x = next(t for t in range(1, p) if (t * t - a) % p == 0)
    q = p
    while q < m:
        q2 = min(q * q, m)
        x = (x - (x * x - a) * pow(2 * x, -1, q2)) % q2
        q = q2
    return x % m


class Elt:
    """a + b sqrt r + c sqrt p + d sqrt(rp), coefficients rationals (exact); p-adic roots of units
    are integers mod p^N, so coefficients are kept as Fractions with p-power-free denominators
    reduced mod p^N at the end."""
    __slots__ = ("c",)
    def __init__(self, c): self.c = c


def mulE(x, y, r, p):
    a, b, c, d = x; e, f, g, h = y
    # basis 1, s = sqrt r, t = sqrt p, st; s^2 = r, t^2 = p
    return (a * e + r * b * f + p * c * g + r * p * d * h,
            a * f + b * e + p * (c * h + d * g),
            a * g + c * e + r * (b * h + d * f),
            a * h + d * e + b * g + c * f)


def sqrt_elt(g, p, r, m):
    """sqrt of the rational integer g (possibly negative) in L, as 4 coefficients (integers mod m
    for the unit part times sqrt r^e1 sqrt p^e2)."""
    v = 0; gg = g
    while gg % p == 0:
        gg //= p; v += 1
    e2 = v % 2; scale = p ** (v // 2)
    u = gg % m
    if pow(u % p, (p - 1) // 2, p) == 1:
        s = hensel_sqrt(u, p, m); e1 = 0
    else:
        s = hensel_sqrt((u * pow(r, -1, m)) % m, p, m); e1 = 1
    coef = [0, 0, 0, 0]
    coef[e1 + 2 * e2] = (s * scale) % m
    return tuple(coef)


def residues(units, F, p, r, signs):
    m = p ** N
    roots = {}
    for g, sg in zip(list(F.gens) + [-1], signs):
        c = sqrt_elt(g, p, r, m)
        roots[g] = tuple((sg * x) % m for x in c)
    B = []
    for i in range(F.dim):
        v = (1, 0, 0, 0)
        for j, g in enumerate(F.gens):
            if i >> j & 1:
                v = tuple(x % m for x in mulE(v, roots[g], r, p))
        B.append(v)
    I = roots[-1]
    out = []
    for u in units:
        acc = [Fr(0)] * 4
        for part, extra in ((u.x.c, (1, 0, 0, 0)), (u.y.c, I)):
            for i, c in enumerate(part):
                c = Fr(c)
                if c:
                    t = mulE(B[i], extra, r, p)
                    for k in range(4):
                        acc[k] += c * (t[k] % m)
        # integrality: coefficients of 1, s integral; of t, st: valuation >= 0 as well (t has v = 1/2)
        if any(x.denominator % p == 0 for x in acc):
            return None
        red = tuple((x.numerator * pow(x.denominator, -1, p)) % p for x in acc[:2])
        out.append(red)
    return out


def colourable(p, S, k, TL):
    V = [(a, b) for a in range(p) for b in range(p)]
    idx = {v: i for i, v in enumerate(V)}
    var = lambda v, c: v * k + c + 1
    from threading import Timer
    with Solver(name="cadical195") as s:
        for i in range(len(V)):
            s.add_clause([var(i, c) for c in range(k)])
        for (a, b) in V:
            for (x, y) in S:
                i, j = idx[(a, b)], idx[((a + x) % p, (b + y) % p)]
                if i < j:
                    for c in range(k):
                        s.add_clause([-var(i, c), -var(j, c)])
        s.add_clause([var(0, 0)])
        tm = Timer(TL, lambda: s.interrupt()); tm.start()
        res = s.solve_limited(expect_interrupt=True); tm.cancel()
        return res


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    F = Field(tuple(d["field_generators"]))
    mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
    units = [mk(xy) for xy in d["units"]]
    p = int(sys.argv[2]); k = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    TL = float(sys.argv[4]) if len(sys.argv) > 4 else 60
    r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
    print(f"{os.path.basename(sys.argv[1])}: {len(units)} units, p = {p}", flush=True)
    seen = {}
    for signs in itertools.product((1, -1), repeat=len(F.gens) + 1):
        res = residues(units, F, p, r, signs)
        if res is None:
            print(f"  signs {signs}: some unit not integral"); continue
        S = frozenset(res)
        if (0, 0) in S:
            print(f"  signs {signs}: a unit reduces to 0 (the place is not where the units live)"); continue
        if S in seen: continue
        seen[S] = colourable(p, S, k, TL)
        print(f"  signs {signs}: {len(S)} residues in F_{p}^2 minus 0; Cay {k}-colourable: {seen[S]}", flush=True)
