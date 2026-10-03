"""Real multiquadratic fields F = Q(sqrt a1, ..., sqrt ak) with no locally constant 4-colouring at any place and no
unit triangle (notes/local_global.md, section 4): sqrt3 is not in F, and every place of F above 2, 3, 7, 11 and 19
contains i, i.e. -1 lies in <a1, ..., ak> Q_p^*2 for p = 2, 3, 7, 11, 19. (The other anisotropic places lie above
primes p = 3 (mod 4) whose levels have Hoffman ratio below 1/4.)

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


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    sq = [d for d in range(2, N) if sqfree(d) == d]
    print("quadratic:", [d for d in sq if admissible([d])][:12])
    bi = sorted({tuple(group_of([a, b])) for a, b in itertools.combinations([d for d in sq if d < 120], 2)
                 if admissible([a, b])}, key=lambda G: (max(G), G))
    print("biquadratic (square classes):", bi[:12])
