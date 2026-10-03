"""Real multiquadratic fields F = Q(sqrt a1, ..., sqrt ak) with no locally constant 4-colouring at any place and no
unit triangle (notes/local_global.md, section 4): sqrt3 is not in F, and every place of F above 2, 3, 7, 11 and 19
contains i, i.e. -1 lies in <a1, ..., ak> Q_p^*2 for p = 2, 3, 7, 11, 19. (Every other place v lies above a prime p
with F_v containing Q_p: either i is in F_v, or p = 3 (mod 4) with p >= 23 and the measurable chromatic number of
Q_p^2, at least 5 there by papers/padic-planes Table 1, bounds chi_loc(v) from below.)

real_place_example() is the smallest prime d = 7 (mod 8) that is a non-residue modulo every prime p = 3 (mod 4) below
100 other than 71: every finite place of Q(sqrt d) then contains i or lies above 71 or a prime p >= 103, where the
measurable chromatic number of Q_p^2 is at least 7 (papers/two-colour-planes, Section 4).

usage: python3 admissible.py [N]   (quadratic d < N, biquadratic with generators below 120)"""
import itertools
import sys


def sqfree(n):
    out, p, m = 1, 2, n
    while p * p <= m:
        e = 0
        while m % p == 0:
            m //= p; e += 1
        if e % 2:
            out *= p
        p += 1
    return out * m


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def group_of(gens):
    """the nontrivial square classes of <gens>, as squarefree positive integers"""
    G = set()
    for r in range(1, len(gens) + 1):
        for c in itertools.combinations(gens, r):
            prod = 1
            for x in c:
                prod *= x
            G.add(sqfree(prod))
    return sorted(G)


def has_i_at(group, p):
    """-1 is in <group> Q_p^*2: some class is -1 times a p-adic square"""
    for x in group:
        if p == 2:
            if x % 8 == 7:          # x odd and -x = 1 (mod 8)
                return True
        elif x % p and legendre(x, p) == -1:   # a unit non-residue, p = 3 (mod 4)
            return True
    return False


BAD = (2, 3, 7, 11, 19)


def admissible(gens):
    G = group_of(gens)
    if 1 in G or 3 in G:
        return False
    return all(has_i_at(G, p) for p in BAD)


def is_prime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))


def real_place_example(skip=(71,), below=100):
    P = [p for p in range(3, below) if is_prime(p) and p % 4 == 3 and p not in skip]
    d = 7
    while not (is_prime(d) and all(legendre(d, p) == -1 for p in P)):
        d += 8
    return d


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    sq = [d for d in range(2, N) if sqfree(d) == d]
    print("quadratic:", [d for d in sq if admissible([d])][:12])
    bi = sorted({tuple(group_of([a, b])) for a, b in itertools.combinations([d for d in sq if d < 120], 2)
                 if admissible([a, b])}, key=lambda G: (max(G), G))
    print("biquadratic (square classes):", bi[:12])
    print("real place: d =", real_place_example(), "(with 71 required too:", real_place_example(skip=()), ")")
