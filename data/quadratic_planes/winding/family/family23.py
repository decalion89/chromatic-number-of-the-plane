"""family23.py d [k]: the configuration of notes/four_colours_23_mod_24.md for d = 23 (mod 24).

U = G u G*u0 u G*conj(u0), where G is the set of the 4(2k+1) rational unit vectors whose denominator divides 5^k and
u0 = ((d-1)/2 + i sqrt d)/((d+1)/2), the unit vector Fischer (1990, Thm 10(ii)) used for additive colourings.
By default k is the least integer >= 1 with d < 21 * 5^(k-1), the bound under which notes/four_colours_11_mod_12.md
(Theorem 1a) proves that no character of ZU maps U into [1/3, 2/3]; the bound is sharp.  Writes in_fam_{d}_{L}.json in the format of certify_w2.py / check_w.py: every vector
((a + b sqrt d) + i (c + e sqrt d))/L as the integer list [a, b, c, e], one per +- pair."""
import sys, json, math
from fractions import Fraction as Fr


def rational_units(N):
    out = set()
    for c in range(1, N + 1):
        if N % c: continue
        for a in range(-c, c + 1):
            b2 = c * c - a * a; b = math.isqrt(b2)
            if b * b == b2:
                for bb in {b, -b}:
                    if math.gcd(math.gcd(abs(a), abs(bb)), c) == 1:
                        out.add((Fr(a, c), Fr(bb, c)))
    return sorted(out)


def least_k(d):
    k = 1
    while d >= 21 * 5 ** (k - 1): k += 1
    return k


def config(d, k):
    assert d % 24 == 23
    m, n, t = (d - 1) // 2, 1, (d + 1) // 2
    assert m * m + d * n * n == t * t
    vecs = []
    for g1, g2 in rational_units(5 ** k):
        vecs.append((g1, Fr(0), g2, Fr(0)))                         # gamma
        for nn in (n, -n):                                           # gamma (m + i nn sqrt d)/t
            vecs.append((g1 * m / t, -g2 * nn / t, g2 * m / t, g1 * nn / t))
    L = 1
    for v in vecs:
        for x in v: L = L * x.denominator // math.gcd(L, x.denominator)
    seen, U = set(), []
    for v in vecs:
        w = tuple(int(x * L) for x in v)
        key = max(w, tuple(-x for x in w))
        if key not in seen:
            seen.add(key); U.append(list(key))
    for a, b, c, e in U:
        assert a * a + d * b * b + c * c + d * e * e == L * L and a * b + c * e == 0
    return L, U


if __name__ == "__main__":
    d = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else least_k(d)
    L, U = config(d, k)
    proven = k >= 1 and d < 21 * 5 ** (k - 1)
    path = f"in_fam_{d}_{L}.json"
    json.dump({"d": d, "D": L, "units": U}, open(path, "w"))
    print(f"d={d} k={k} N={5**k}: L={L}, {len(U)} vectors up to sign -> {path}; "
          f"the note's bound d < 21 * 5^(k-1) {'holds' if proven else 'FAILS'}")
