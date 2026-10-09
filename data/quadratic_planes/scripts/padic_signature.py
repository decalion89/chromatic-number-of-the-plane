"""Does a colouring of a growth graph over Q(sqrt d) look like a reduction at a place above p?

For each place v above p where d is a square modulo p (two places, one for each square root of d), and each
level r, the points of the graph are grouped by their residue modulo p^r at v, and the script counts, over all
pairs of points in the same class, how often the two colours agree. A colouring that factors through level r at
v gives agreement 1 there; a colouring unrelated to v gives about the sum of the squared colour frequencies.

Caution: colours agree more often at even distance in the graph, in any colouring. If the directions collide
modulo p^r, many pairs of one class are close in the graph and the agreement rises for that reason alone, as it
did at 19 for Q(sqrt47) with D = 240 (docs/research-log.md, 3 October). --by-distance splits the pairs by their
distance in the graph (1, 2, 3, 4, 5 or more); only the long-range agreement says something about the place.

usage: python3 padic_signature.py STATE.json[.gz] D_SQUAREFREE p1,p2,... [levels] [--by-distance]
       (STATE has "P": points (a, b, c, e) meaning ((a + b sqrt d)/D, (c + e sqrt d)/D), "U": the directions,
        "col", and "D")"""
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


def by_distance(P, U, col, D, d, p, sign, r):
    """pairs of one residue class, split by graph distance (5 means 5 or more): {k: (pairs, agreements)}"""
    idx = {tuple(q): i for i, q in enumerate(P)}
    nb = []
    for q in P:
        nb.append([idx[t] for t in (tuple(a + b for a, b in zip(q, u)) for u in U) if t in idx])
    ball2 = {}

    def b2(i):
        if i not in ball2:
            s = set(nb[i]); s.add(i)
            for j in nb[i]:
                s.update(nb[j])
            ball2[i] = s
        return ball2[i]

    mod = p ** r
    s = sign * sqrt_mod_pk(d, p, r) % mod
    Dinv = pow(D, -1, mod)
    cls = defaultdict(list)
    for i, (a, b, c, e) in enumerate(P):
        cls[((a + b * s) * Dinv % mod, (c + e * s) * Dinv % mod)].append(i)
    out = defaultdict(lambda: [0, 0])
    for v in cls.values():
        for x in range(len(v)):
            for y in range(x + 1, len(v)):
                i, j = v[x], v[y]
                ni, nj = set(nb[i]), set(nb[j])
                if j in ni:
                    k = 1
                elif ni & nj:
                    k = 2
                elif b2(i) & nj:
                    k = 3
                elif b2(i) & b2(j):
                    k = 4
                else:
                    k = 5
                out[k][0] += 1
                out[k][1] += col[i] == col[j]
    return dict(out)


if __name__ == "__main__":
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path, d = args[0], int(args[1])
    primes = [int(x) for x in args[2].split(",")]
    levels = [int(x) for x in args[3].split(",")] if len(args) > 3 else [1, 2, 3]
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
                if "--by-distance" in flags:
                    for k, (n, a) in sorted(by_distance(P, st["U"], col, D, d, p, sign, r).items()):
                        print(f"    distance {'5 or more' if k == 5 else k}: {n} pairs, agreement {a / n:.3f}",
                              flush=True)
