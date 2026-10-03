"""Does a colouring of a growth graph over Q(sqrt d) look like a reduction at a place above p?

For each place v above p where d is a square modulo p (two places, one for each square root of d), and each
level r, the points of the graph are grouped by their residue modulo p^r at v, and the script counts, over all
pairs of points in the same class, how often the two colours agree. A colouring that factors through level r at
v gives agreement 1 there; a colouring unrelated to v gives about the sum of the squared colour frequencies.

usage: python3 padic_signature.py STATE.json[.gz] D_SQUAREFREE p1,p2,... [levels]
       (STATE has "P": points (a, b, c, e) meaning ((a + b sqrt d)/D, (c + e sqrt d)/D), "col", and "D")"""
import gzip
import json
import sys
from collections import defaultdict


def sqrt_mod_pk(a, p, k):
    """a square root of a modulo p^k by Hensel lifting (p odd, a a nonzero square modulo p)"""
    r = next(x for x in range(1, p) if (x * x - a) % p == 0)
    mod = p
    for _ in range(1, k):
        mod *= p
        r = (r - (r * r - a) * pow(2 * r, -1, mod)) % mod
    return r


def signature(P, col, D, d, p, sign, r):
    mod = p ** r
    s = sign * sqrt_mod_pk(d, p, r) % mod
    Dinv = pow(D, -1, mod)
    cls = defaultdict(lambda: [0] * (max(col) + 1))
    for (a, b, c, e), k in zip(P, col):
        cls[((a + b * s) * Dinv % mod, (c + e * s) * Dinv % mod)][k] += 1
    pairs = agree = 0
    for h in cls.values():
        n = sum(h)
        pairs += n * (n - 1) // 2
        agree += sum(x * (x - 1) // 2 for x in h)
    return len(cls), pairs, agree


if __name__ == "__main__":
    path, d = sys.argv[1], int(sys.argv[2])
    primes = [int(x) for x in sys.argv[3].split(",")]
    levels = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [1, 2, 3]
    st = json.load(gzip.open(path) if path.endswith(".gz") else open(path))
    P, col, D = st["P"], st["col"], st.get("D", 240)
    freq = [col.count(k) / len(col) for k in range(max(col) + 1)]
    print(f"{len(P)} points, D = {D}; agreement of a colouring unrelated to any place: {sum(f * f for f in freq):.3f}")
    for p in primes:
        if pow(d % p, (p - 1) // 2, p) != 1:
            print(f"p = {p}: d is not a nonzero square modulo p, so no place above p has F_v = Q_p; skipped")
            continue
        for sign in (1, -1):
            for r in levels:
                ncls, pairs, agree = signature(P, col, D, d, p, sign, r)
                print(f"p = {p}, place {'+' if sign > 0 else '-'}, level {r}: {ncls} classes, {pairs} pairs in a "
                      f"class, agreement {agree / pairs if pairs else float('nan'):.3f}", flush=True)
