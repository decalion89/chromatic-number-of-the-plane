"""Is rigidity emerging as the grown graphs get bigger?

Every coset colouring colours x + u and x + v alike whenever u - v lies in 5M (kappa-pairs are
the prime example: kappa*u - u is in 5M).  Take a coloured checkpoint and, over every vertex x
and every two of its unit neighbours x + u, x + v, compare the alike-rate of the 5M pairs with
the alike-rate of all other non-adjacent pairs.  In a colouring that is coset on a region the
first rate is 1; in a colouring that does not feel the rigidity the two rates are about equal.

usage: python3 kpairs.py <checkpoint.json> [more checkpoints ...]"""
import sys, json, os
from fractions import Fraction as Fr
from math import gcd
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])


def analyse(path):
    d = json.load(open(path)); col = d.get("colouring")
    assert col, "checkpoint has no colouring"
    F = Field(tuple(d["field_generators"]))
    mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
    V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]; A = d["A"]
    Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
    raw = [tuple((p - V[A]).x.c) + tuple((p - V[A]).y.c) for p in V]
    den = 1
    for v in Ev + raw:
        for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev)})
    B = echelon(E)
    iv = lambda v: tuple(int(t) for t in coords(B, tuple(int(Fr(x) * den) for x in v)))
    cu = [iv(v) for v in Ev]; cp = [iv(v) for v in raw]
    idx = {c: i for i, c in enumerate(cp)}
    nu = len(U); uset = set(cu)
    # classify unit pairs: adjacent (difference is a unit), 5M, other; bucket by |u - v|^2
    cls = {}
    for i in range(nu):
        for j in range(i + 1, nu):
            dv = tuple(a - b for a, b in zip(cu[i], cu[j]))
            if dv in uset or not any(dv): continue
            w = U[i] - U[j]; d2 = round(w.fx ** 2 + w.fy ** 2, 6)
            cls[(i, j)] = (all(t % 5 == 0 for t in dv), d2)
    tally = defaultdict(lambda: [0, 0])      # (in5, d2) -> [alike, total]
    for x, c in enumerate(cp):
        nb = [i for i in range(nu) if tuple(a + b for a, b in zip(c, cu[i])) in idx]
        ci = [col[idx[tuple(a + b for a, b in zip(c, cu[i]))]] for i in nb]
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                k = cls.get((nb[a], nb[b]))
                if k is None: continue
                t = tally[k]; t[0] += ci[a] == ci[b]; t[1] += 1
    tot = {True: [0, 0], False: [0, 0]}
    for (in5, d2), (al, n) in tally.items(): tot[in5][0] += al; tot[in5][1] += n
    print(f"{os.path.basename(path)}: n={len(V)}, {nu} units, rank {len(B)}")
    for in5 in (True, False):
        al, n = tot[in5]
        print(f"  {'5M pairs   ' if in5 else 'other pairs'}: {al}/{n} alike = {al / max(n, 1):.3f}")
    # the 5M distance classes one by one, against the other pairs at the same distance
    for d2 in sorted({k[1] for k in tally if k[0]}):
        a5, n5 = tally[(True, d2)]; ao, no = tally.get((False, d2), [0, 0])
        print(f"    |u-v|^2 = {d2:<9}: 5M {a5}/{n5} = {a5 / max(n5, 1):.3f};  other {ao}/{no} = {ao / max(no, 1):.3f}")


for p in sys.argv[1:]: analyse(p)
