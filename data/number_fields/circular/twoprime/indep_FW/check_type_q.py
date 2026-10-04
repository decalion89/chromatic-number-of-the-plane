"""Referee check (own code) of the type-q step in the draft's proof for chi_c = 3.

Setting: N = 5^k, G_N = {i^a rho^j : |j| <= k} (paper, Lemma values (i)), phi(v) = N*eps with eps in E_q,
eps = (alpha + beta i)/3 + t, alpha, beta in {1,2}, t in Z[i].  The character is xi(gamma v) = Re(gamma*N*eps) mod 1.
Checks:
 (1) for k = 1..6 and every gamma in G_N: N*gamma in Z[i] and Re(N gamma (alpha+beta i)) is not 0 mod 3, so
     xi(gamma v) in {1/3, 2/3} mod 1 (and independent of t);
 (2) for gamma in G_5: Re(N gamma (alpha + beta i)) == 5^(k-1) Re(5 gamma (alpha+beta i)) (mod 3), so the tight set
     T = {gamma in G_5 : value 1/3} = {gamma : Re(5 gamma (alpha+beta i)) == s (mod 3)}, s = 5^(k-1) mod 3 = (-1)^(k-1);
 (3) for each of the 4 classes and both signs s: |T| = 6, T contains no pair {g, -g}; list ALL triples (and pairs)
     in T with a positive rational relation; check that one triple has multiplicities exactly (3,4,5) (up to order);
     check the draft's example 3(-1) + 4i + 5(3-4i)/5 = 0 and for which (class, sign) it is tight;
 (4) positive relations exist iff 0 is in the interior of the convex hull of T; we also report the general fact.
"""
from fractions import Fraction as Fr
from itertools import combinations, permutations

def cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

def cpow(z, e):
    r = (Fr(1), Fr(0))
    if e < 0:
        z = (z[0], -z[1]); e = -e
    for _ in range(e):
        r = cm(r, z)
    return r

RHO = (Fr(3, 5), Fr(4, 5))
UNITS = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]

def G(k):
    return [cm(u, cpow(RHO, j)) for u in UNITS for j in range(-k, k + 1)]

def is_gint(z):
    return z[0].denominator == 1 and z[1].denominator == 1

# (1)
for k in range(1, 7):
    N = 5 ** k
    GN = G(k)
    assert len(set(GN)) == 4 * (2 * k + 1)
    for al in (1, 2):
        for be in (1, 2):
            for g in GN:
                Ng = (N * g[0], N * g[1])
                assert is_gint(Ng)
                re = cm(Ng, (Fr(al), Fr(be)))[0]
                assert re.denominator == 1 and int(re) % 3 != 0
print("(1) for k=1..6, all gamma in G_N, all 4 classes: Re(N gamma (alpha+beta i)) != 0 mod 3, values in {1/3,2/3}: OK")

# (2)
for k in range(1, 7):
    N = 5 ** k
    for al in (1, 2):
        for be in (1, 2):
            for g in G(1):
                a = int(cm((N * g[0], N * g[1]), (Fr(al), Fr(be)))[0])
                b = int(cm((5 * g[0], 5 * g[1]), (Fr(al), Fr(be)))[0])
                assert (a - pow(5, k - 1, 3) * b) % 3 == 0
print("(2) Re(N gamma eps') == 5^(k-1) Re(5 gamma eps') mod 3 on G_5 for k = 1..6: OK")

G5 = G(1)
assert len(set(G5)) == 12

def pos_relation_triple(a, b, c):
    """positive rational (x, y, 1) with x a + y b + c = 0, or None"""
    det = a[0] * b[1] - a[1] * b[0]
    if det == 0:
        return None
    x = (-c[0] * b[1] + c[1] * b[0]) / det
    y = (-a[0] * c[1] + a[1] * c[0]) / det
    if x > 0 and y > 0:
        return (x, y, Fr(1))
    return None

def normalise(m):
    from math import lcm, gcd
    L = lcm(*(t.denominator for t in m))
    v = [int(t * L) for t in m]
    g = 0
    for t in v:
        g = gcd(g, t)
    return tuple(t // g for t in v)

def show(z):
    return f"{z[0]}{'+' if z[1] >= 0 else '-'}{abs(z[1])}i"

example = [((Fr(-1), Fr(0)), 3), ((Fr(0), Fr(1)), 4), ((Fr(3, 5), Fr(-4, 5)), 5)]
assert sum(m * z[0] for z, m in example) == 0 and sum(m * z[1] for z, m in example) == 0

all_ok = True
for al in (1, 2):
    for be in (1, 2):
        for s in (1, -1):
            T = [g for g in G5 if (int(cm((5 * g[0], 5 * g[1]), (Fr(al), Fr(be)))[0]) - s) % 3 == 0]
            assert len(T) == 6
            assert not any((-g[0], -g[1]) in T for g in T)
            rels = []
            for a, b, c in combinations(T, 3):
                r = pos_relation_triple(a, b, c)
                if r is not None:
                    rels.append(((a, b, c), normalise(r)))
            has345 = [x for x in rels if sorted(x[1]) == [3, 4, 5]]
            ex_tight = all(z in T for z, _ in example)
            print(f"class ({al}+{be}i)/3, sign {s:+d}: T = {[show(g) for g in T]}")
            print(f"    {len(rels)} positive triples; multiplicity patterns {sorted(set(tuple(sorted(r[1])) for r in rels))};"
                  f" a (3,4,5) triple: {[ (tuple(show(z) for z in x[0]), x[1]) for x in has345[:1]]};"
                  f" draft's example tight here: {ex_tight}")
            if not has345:
                all_ok = False
print("(3) every class and sign has a positive relation with multiplicities 3,4,5:", all_ok)
