"""The minimum rank at which unit directions can block a coset 5-colouring.

Rank 2 cannot: PG(1,5)'s hyperplanes are single points, so covering needs six
directions, and a plane lattice's unit circle holds at most six points hence
three.  Rank 6 can: Q(zeta_7)'s modulus-one elements do it with 87.  The gap
in between is worth closing, because rank is what costs density -- a discrete
module has rank at most 2, and every step up makes the ball's boundary
dominate its interior further.

A pencil, the cheapest cover, is out of reach at every rank: six concurrent
hyperplanes means the six directions lie in one 2-dimensional subspace, and
there they must hit all six points of a PG(1,5) while sharing one norm, which
no square class allows.  So any blocking set of unit vectors is NOT minimal,
and the question is how many directions a given rank needs.

This searches subsets of Q(zeta_7)'s modulus-one steps by the rank of the
module they generate.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism

F = CycloField(7)
D = F.degree
one = F.rational(1)


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def rank(vs):
    M = [list(map(Fraction, v)) for v in vs]
    if not M:
        return 0
    r = 0
    for c in range(len(M[0])):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r


steps = set()
for coeffs in itertools.product(range(-3, 4), repeat=4):
    a = tuple(Fraction(c) for c in coeffs) + (Fraction(0), Fraction(0))
    if not any(a):
        continue
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) != one:
        continue
    d = 1
    for x in u:
        d = d * x.denominator // gcd(d, x.denominator)
    if d % 5:
        steps.add(u)
steps = sorted(steps)
print(f"{len(steps)} modulus-one steps", flush=True)


def directions(sub):
    den = 1
    for v in sub:
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
    out = set()
    for v in sub:
        w = tuple(int(x * den) for x in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        w = tuple(t // g for t in w) if g > 1 else w
        out.add(min(w, tuple(-t for t in w)))
    return sorted(out)


# Random subsets of Q(zeta_7)'s steps have full rank almost always, so the
# low-rank cases have to be CONSTRUCTED: span a submodule with a few steps,
# then collect every unit step that lands inside it.  That is the whole unit
# circle of that submodule, which is what a graph built there would have.
rng = random.Random(17)
best = {}
for target in (3, 4, 5):
    seen = 0
    for trial in range(300):
        seed = rng.sample(steps, target)
        if rank(seed) != target:
            continue
        inside = [u for u in steps if rank(seed + [u]) == target]
        if len(inside) < 6:
            continue
        seen += 1
        dirs = directions(inside)
        if rank(dirs) != target:
            continue
        blocks = has_homomorphism(dirs, 5)[0] is None
        cur = best.get(target)
        if blocks and (cur is None or not cur[0]):
            best[target] = (True, len(dirs), len(inside))
            print(f"  rank {target}: BLOCKS -- {len(dirs)} directions from "
                  f"{len(inside)} unit steps", flush=True)
            break
        if cur is None or len(dirs) > cur[1]:
            best[target] = (False, len(dirs), len(inside))
    if target not in best:
        print(f"  rank {target}: no submodule with 6+ unit steps found",
              flush=True)
    elif not best[target][0]:
        b, nd, nu = best[target]
        print(f"  rank {target}: no block in {seen} submodules; richest had "
              f"{nd} directions from {nu} unit steps", flush=True)

print("\nsmallest blocking rank found:",
      min([r for r, v in best.items() if v[0]], default="none at rank <= 5"))
print("  (rank 2 is impossible by the covering count; rank 6 blocks with 87)")
