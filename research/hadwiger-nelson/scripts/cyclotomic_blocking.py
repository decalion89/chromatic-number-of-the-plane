"""Make 5 inert: the multiquadratic families were the wrong shape.

The projective direction count of the modulus-one group mod 5 is a product of
small factors -- 3 where 5 is inert in a quadratic factor, 2 where it splits --
so a multiquadratic field of degree 2t gives at most 3^t directions against
PG(2t-1,5)'s (5^2t - 1)/4 points. Measured at 27 against 97656. The ratio
(3/25)^t collapses, and no amount of extra generators helps: the group is what
it is.

The cause is the splitting. In Q(sqrt-d1,...,sqrt-dt) the Galois group is
(Z/2)^t and the decomposition group at 5 is cyclic, so the residue degree is
at most 2 and O/5 breaks into many tiny fields. To get a LARGE norm-one group
the decomposition group must be everything, so the Galois group has to be
cyclic and 5 inert -- which is Q(zeta_n) with 5 a primitive root mod n.

n = 7 qualifies: 5 has order 6 in (Z/7)^*. Then O/5 = F_5^6 and the norm-one
subgroup has (5^6 - 1)/(5^3 - 1) = 126 elements against PG(5,5)'s 3906 points.

Modulus-one elements come free from Hilbert 90: u = alpha / conj(alpha) has
modulus one for every alpha, and its denominator is a norm, which only has to
be kept coprime to 5.
"""
import sys, time, random, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from math import gcd
from hn.homcol import has_homomorphism

N = 7
D = N - 1          # degree of Q(zeta_7)


def mul(x, y):
    """Multiply in Z[x]/(1 + x + ... + x^6), i.e. Z[zeta_7]."""
    raw = [0] * (2 * D + 1)
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                if b:
                    raw[i + j] += a * b
    # reduce x^k for k >= D using x^6 = -(1 + x + ... + x^5)
    for k in range(2 * D, D - 1, -1):
        c = raw[k]
        if c:
            raw[k] = 0
            for j in range(D):
                raw[k - D + j] -= c
    return tuple(raw[:D])


def conj(x):
    """zeta -> zeta^-1 = zeta^6, then reduce."""
    out = [0] * D
    for i, a in enumerate(x):
        if not a:
            continue
        k = (-i) % N
        if k < D:
            out[k] += a
        else:                       # zeta^6 = -(1 + ... + zeta^5)
            for j in range(D):
                out[j] -= a
    return tuple(out)


def norm_to_int(x):
    """x * conj(x) is real; return it as a rational-free integer if it is."""
    p = mul(x, conj(x))
    return p


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


def norm_dir(v, p=5):
    v = tuple(t % p for t in v)
    if not any(v):
        return None
    for lam in range(1, p):
        w = tuple((t * lam) % p for t in v)
        for i in range(len(w)):
            if w[i] == 1 and all(w[j] == 0 for j in range(i)):
                return w
    return None


# u = alpha / conj(alpha): represent as the pair, clear by multiplying by
# conj(alpha)^-1 mod 5 instead of exactly -- all we need is the reduction.
def inv_mod5(x):
    """Inverse in (Z[zeta]/5)^*, by brute-force extended power."""
    # the group has order 5^6 - 1; x^(5^6 - 2) is the inverse
    e = 5 ** D - 2
    r, b = tuple([1] + [0] * (D - 1)), tuple(t % 5 for t in x)
    while e:
        if e & 1:
            r = tuple(t % 5 for t in mul(r, b))
        b = tuple(t % 5 for t in mul(b, b))
        e >>= 1
    return r


random.seed(7)
t0 = time.time()
dirs, vecs = set(), []
tried = 0
for _ in range(4000):
    a = tuple(random.randint(-3, 3) for _ in range(D))
    if not any(a):
        continue
    ca = conj(a)
    if all(t % 5 == 0 for t in ca):
        continue
    tried += 1
    u = tuple(t % 5 for t in mul(a, inv_mod5(ca)))
    d = norm_dir(u)
    if d and d not in dirs:
        dirs.add(d)
        vecs.append(u)
print(f"Q(zeta_7): {tried} alphas -> {len(dirs)} distinct directions mod 5, "
      f"rank {rank_mod(vecs, 5)} of {D}  [{time.time()-t0:.0f}s]", flush=True)
print(f"  PG({D-1},5) has {(5**D - 1)//4} points; the norm-one subgroup has "
      f"{(5**D - 1)//(5**(D//2) - 1)} elements", flush=True)
phi, why = has_homomorphism(vecs, 5)
ok = phi is not None and all(
    sum(p * t for p, t in zip(phi, v)) % 5 for v in vecs)
print(f"  Z/5 homomorphism: " + (f"exists (verified {ok})" if phi
                                 else "*** NONE -- BLOCKED ***"), flush=True)
