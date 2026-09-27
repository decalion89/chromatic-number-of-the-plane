"""Residue gate: colourings of a module through one place above an unramified prime p.

The unit vectors u = x + iy of the module are embedded in Q_{p^2} = Q_p(sqrt r) (every square root
of the field's generators and of -1 exists there), for every choice of signs, i.e. at every place
of K above p.  If all unit vectors are p-integral there, the residues S = {u mod p} lie in F_{p^2}
and any colouring of Cay(F_{p^2}, S) pulls back to a colouring of every graph along these units
(z -> residue of z relative to its coset of the integers, as in hn/adelic.py).  So:
  chi(Cay(F_{p^2}, S)) <= 5  ==>  the module is 5-colourable, whatever graph is grown in it.
Units with p in the denominator make the place non-integral and are reported instead.

usage: python3 residuegate.py <module.json> <p1,p2,...> [k=5]"""
import sys, os, json, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from hn.field import Field
from hn.geometry import Point
from pysat.solvers import Solver

PREC = 12
TL = float(os.environ.get('TL', '60'))   # seconds per colourability test


def sqrt_in_Qp2(g, p, r, N=PREC):
    """A square root of the rational integer g in Q_p(sqrt r) as (a, b) mod p^N (value a + b sqrt r),
    for p not dividing g, p odd, r a non-residue."""
    m = p ** N
    if pow(g % p, (p - 1) // 2, p) == 1:          # g is a square in Q_p
        s = next(x for x in range(1, p) if (x * x - g) % p == 0)
        q = p
        while q < m:
            q2 = q * q if q * q <= m else m
            s = (s - (s * s - g) * pow(2 * s, -1, q2)) % q2
            q = q2
        return (s % m, 0)
    # g / r is a square: sqrt g = t sqrt r with t^2 = g/r
    h = (g * pow(r, -1, m)) % m
    t = next(x for x in range(1, p) if (x * x - h) % p == 0)
    q = p
    while q < m:
        q2 = q * q if q * q <= m else m
        t = (t - (t * t - h) * pow(2 * t, -1, q2)) % q2
        q = q2
    return (0, t % m)


def mul(u, v, r, m):
    return ((u[0] * v[0] + r * u[1] * v[1]) % m, (u[0] * v[1] + u[1] * v[0]) % m)


def residues(units, F, p, r, signs, N=PREC):
    """Residues mod p of the units at the place given by the sign choices; None if some unit is
    not p-integral there."""
    m = p ** N
    roots = {}
    for g, sg in zip(list(F.gens) + [-1], signs):
        a, b = sqrt_in_Qp2(g % m if g > 0 else (m + g), p, r, N)
        roots[g] = ((sg * a) % m, (sg * b) % m)
    def basis_elt(prod_mask):
        v = (1, 0)
        for i, g in enumerate(F.gens):
            if prod_mask >> i & 1:
                v = mul(v, roots[g], r, m)
        return v
    B = [basis_elt(i) for i in range(F.dim)]
    I = roots[-1]
    out = []
    for u in units:
        acc = (0, 0)
        for part, extra in ((u.x.c, (1, 0)), (u.y.c, I)):
            for i, c in enumerate(part):
                c = Fr(c)
                if c == 0:
                    continue
                if c.denominator % p == 0:
                    return None
                cc = (c.numerator * pow(c.denominator, -1, m)) % m
                t = mul(B[i], extra, r, m)
                acc = ((acc[0] + cc * t[0]) % m, (acc[1] + cc * t[1]) % m)
        out.append((acc[0] % p, acc[1] % p))
    return out


def colourable(p, S, k):
    V = [(a, b) for a in range(p) for b in range(p)]
    idx = {v: i for i, v in enumerate(V)}
    var = lambda v, c: v * k + c + 1
    with Solver(name="cadical195") as s:
        for i in range(len(V)):
            s.add_clause([var(i, c) for c in range(k)])
        for (a, b) in V:
            for (x, y) in S:
                j, i = idx[((a + x) % p, (b + y) % p)], idx[(a, b)]
                if i < j:
                    for c in range(k):
                        s.add_clause([-var(i, c), -var(j, c)])
        s.add_clause([var(0, 0)])
        # hard instances (the finite planes at non-split places) are cut off: None = unknown
        from threading import Timer
        timer = Timer(TL, lambda: s.interrupt()); timer.start()
        r = s.solve_limited(expect_interrupt=True)
        timer.cancel()
        return r


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    F = Field(tuple(d["field_generators"]))
    mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
    if "units" in d and len(d["units"]) > 1:
        units = [mk(xy) for xy in d["units"]]
    else:
        from hn.graph import build_graph
        P = [mk(xy) for xy in d["points"]]
        g = build_graph(P)
        U = {}
        for a, b in g.edges():
            for w in (P[b] - P[a], P[a] - P[b]):
                U.setdefault(w, w)
        units = list(U)
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    print(f"{os.path.basename(sys.argv[1])}: {len(units)} unit vectors, field {F.gens}", flush=True)
    for p in [int(x) for x in sys.argv[2].split(",")]:
        if any(g % p == 0 for g in F.gens) or p == 2:
            print(f"  p = {p}: ramified or 2, skipped"); continue
        r = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        seen = {}
        for signs in itertools.product((1, -1), repeat=len(F.gens) + 1):
            res = residues(units, F, p, r, signs)
            if res is None:
                seen.setdefault("non-integral", 0); seen["non-integral"] += 1
                continue
            S = frozenset(res)
            if (0, 0) in S:
                seen.setdefault("zero residue", 0); seen["zero residue"] += 1
                continue
            if S not in seen:
                seen[S] = colourable(p, S, k)
                print(f"    place {signs}: {len(S)} residues, {k}-colourable: {seen[S]}", flush=True)
        summ = {("non-integral" if s == "non-integral" else "zero residue" if s == "zero residue" else f"|S|={len(s)}"): v for s, v in seen.items()}
        bad = [s for s, v in seen.items() if not isinstance(s, str) and v]
        print(f"  p = {p}: " + ", ".join(f"{a}: {b}" for a, b in summ.items())
              + (f"   ==> {k}-COLOURABLE through a place above {p}" if bad else ""), flush=True)
