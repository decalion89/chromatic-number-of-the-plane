"""Do the directions form a projective line? That is the whole question.

The unit vectors of Q(sqrt-d1, sqrt-d2) are the products of modulus-one
elements of each factor, and their images mod 5 saturate: u1 and u2 have
finite order there, so no matter how many products are taken only finitely
many DIRECTIONS appear -- six, for (-7, -11), however far the exponents run.

Six is exactly the size of the smallest blocking set in PG(3,5): a projective
line. So if the six directions were the six points of a line, every phi would
kill one and no coset colouring could exist. They are not, or a phi would not
have been found -- so the question is which pairs of fields give directions
that DO span only rank two and cover their line.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from math import gcd, isqrt
from hn.homcol import has_homomorphism


def units_of(d):
    """Modulus-one generators (p + q sqrt(-d))/m with p^2 + d q^2 = m^2."""
    out = []
    for m in range(2, 40):
        for q in range(1, m + 1):
            r = m * m - d * q * q
            if r <= 0:
                continue
            p = isqrt(r)
            if p * p == r and gcd(gcd(p, q), m) == 1 and p and q and m % 5:
                # A denominator divisible by 5 breaks the reduction: the
                # module stops being 5-integral, clearing denominators makes
                # every vector divisible by 5, and the search reports "no
                # homomorphism" for a reason that is pure arithmetic bookkeeping.
                # That false positive is exactly what this filter prevents.
                out.append((Fr(p, m), Fr(q, m)))
    return out


def field_mul(A, B, x, y):
    return (x[0]*y[0] + A*x[1]*y[1] + B*x[2]*y[2] + A*B*x[3]*y[3],
            x[0]*y[1] + x[1]*y[0] + B*(x[2]*y[3] + x[3]*y[2]),
            x[0]*y[2] + x[2]*y[0] + A*(x[1]*y[3] + x[3]*y[1]),
            x[0]*y[3] + x[3]*y[0] + x[1]*y[2] + x[2]*y[1])


def norm_dir(v, p=5):
    v = tuple(t % p for t in v)
    if not any(v):
        return None
    for lam in range(1, p):
        w = tuple((t * lam) % p for t in v)
        for i in range(len(w)):
            if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                return w
    return v


def rank_mod(vecs, p=5):
    rows = [list(v) for v in vecs]
    if not rows:
        return 0
    n, r = len(rows[0]), 0
    for c in range(n):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] % p), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c] % p, p - 2, p)
        rows[r] = [(v * inv) % p for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] % p:
                f = rows[i][c] % p
                rows[i] = [(a - f * b) % p for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def integerise(vecs):
    den = 1
    for v in vecs:
        for q in v:
            den = den * q.denominator // gcd(den, q.denominator)
    ints = {tuple(int(q * den) for q in v) for v in vecs}
    c = 0
    for v in ints:
        for t in v:
            c = gcd(c, abs(t))
    return sorted({tuple(t // c for t in v) for v in ints}) if c > 1 else sorted(ints)


t0 = time.time()
DS = [7, 11, 15, 19, 23, 31, 39, 43, 47, 55, 59, 67, 71, 79, 83, 87, 91, 95]
print(f"{'d1':>4} {'d2':>4} {'units':>6} {'dirs':>5} {'rank':>5}  Z/5",
      flush=True)
best = None
for d1, d2 in itertools.combinations(DS, 2):
    g1, g2 = units_of(d1)[:3], units_of(d2)[:3]
    if not g1 or not g2:
        continue
    A, B = -d1, -d2
    one = (Fr(1), Fr(0), Fr(0), Fr(0))
    elems = {one}
    gens = [(p, q, Fr(0), Fr(0)) for p, q in g1] + \
           [(p, Fr(0), q, Fr(0)) for p, q in g2]
    for _ in range(3):
        new = set(elems)
        for e in elems:
            for g in gens:
                new.add(field_mul(A, B, e, g))
        elems = new
        if len(elems) > 900:
            break
    iv = integerise(sorted(elems))
    dirs = {norm_dir(v) for v in iv} - {None}
    r = rank_mod(iv, 5)
    phi, why = has_homomorphism(iv, 5)
    zero = [v for v in iv if all(t % 5 == 0 for t in v)]
    flag = ("   [DEGENERATE: %d vectors are 0 mod 5]" % len(zero)) if zero \
        else ("" if phi else "   *** NO HOMOMORPHISM ***")
    if not phi or (best is None or len(dirs) > best[0]):
        best = (len(dirs), d1, d2)
        print(f"{d1:>4} {d2:>4} {len(iv):>6} {len(dirs):>5} {r:>5}  "
              + ("exists" if phi else "NONE") + flag
              + f"   [{time.time()-t0:.0f}s]", flush=True)
    if not phi:
        break
print(f"best direction count: {best}  [{time.time()-t0:.0f}s]", flush=True)
