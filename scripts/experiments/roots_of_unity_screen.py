"""Roots of unity as unit steps: the best ratio of steps to module rank.

A coset colouring phi: M -> Z/5 exists unless the hyperplanes d-perp, one per
edge vector, cover the whole dual (Z/5)^r. A counting estimate says a random
phi survives m hyperplanes with probability (4/5)^m, so coverage needs roughly

    m  >  r * log 5 / log(5/4)  =  7.2 r

distinct unit vectors per unit of module rank. The three-hexagon gadget has
r = 4 and m = 18, needing 29 -- and indeed a phi exists. de Grey's G has
r = 16 and m = 133 against a threshold of 115, marginal, and a phi exists too.

Z[zeta_n] maximises the ratio for free: 2n unit vectors -- every root of
unity is a unit vector of the plane -- against rank phi(n). The ratio
2n/phi(n) is large exactly when n has many small prime factors, reaching
8.75 at n = 210. Whether that actually blocks is a SAT question.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import has_homomorphism


def polydiv(a, b):
    """Exact division of integer polynomials, low-degree-first coefficients."""
    a = list(a)
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(q) - 1, -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j, bj in enumerate(b):
            a[i + j] -= c * bj
    assert all(v == 0 for v in a), "not an exact division"
    return q


_CYC = {}


def cyclotomic(n):
    """Phi_n, from x^n - 1 = prod over d | n of Phi_d."""
    if n in _CYC:
        return _CYC[n]
    num = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            num = polydiv(num, cyclotomic(d))
    _CYC[n] = num
    return num


def totient(n):
    return len([k for k in range(1, n + 1)
                if all((n % p or k % p) for p in range(2, n + 1)
                       if n % p == 0)]) if False else sum(
        1 for k in range(1, n + 1) if _gcd(k, n) == 1)


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def roots_as_vectors(n):
    """zeta^k for k = 0..n-1 in the power basis of Z[zeta_n]."""
    phi = cyclotomic(n)
    d = len(phi) - 1
    out = []
    for k in range(n):
        r = [0] * (k + 1)
        r[k] = 1
        for i in range(len(r) - 1, d - 1, -1):
            c = r[i]
            if c:
                for j, pj in enumerate(phi):
                    r[i - d + j] -= c * pj
        r = r + [0] * (d - len(r))
        out.append(tuple(r[:d]))
    return out, d


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


t0 = time.time()
print(f"{'n':>4} {'phi(n)':>7} {'steps':>6} {'rank5':>6} {'2n/phi':>7}  Z/5", flush=True)
for n in [3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 18, 20, 21, 24, 28, 30, 35,
          36, 40, 42, 60, 70, 84, 105, 120, 140, 210]:
    vecs, d = roots_as_vectors(n)
    # differences of roots are the real edge vectors, but a step SET is what
    # a graph is built from: the steps themselves are the edge vectors.
    uniq = sorted({v for v in vecs if any(v)})
    r5 = rank_mod(uniq, 5)
    phi, why = has_homomorphism(uniq, 5)
    ok = phi is not None and all(
        sum(p * v for p, v in zip(phi, w)) % 5 for w in uniq)
    print(f"{n:>4} {totient(n):>7} {len(uniq):>6} {r5:>6} "
          f"{2*n/totient(n):>7.2f}  "
          + (f"exists (verified {ok})" if phi else f"NONE -- {why}")
          + f"   [{time.time()-t0:.0f}s]", flush=True)
