"""Gamma = Z/8, S = {1,2,6,7}, c(2k) = k, c(2k+1) = k + 2 (mod 4): proper, has the tight 4-cycle
0 -> 2 -> 4 -> 6 -> 0 (steps 2,2,2,2; not a square), no tight square, A = B everywhere,
a(s)/4 is the character x -> 5x/8, kappa(S) = 1/4, and the tight digraph has a cycle."""
from fractions import Fraction as Fr
n, S = 8, [1, 2, 6, 7]
c = {x: (x // 2 + (2 if x % 2 else 0)) % 4 for x in range(n)}
ell = {(x, s): (c[(x + s) % n] - c[x]) % 4 for x in range(n) for s in S}
assert all(v for v in ell.values()), "improper"
sq = [(x, s, t) for x in range(n) for s in S for t in S
      if ell[(x, s)] == 1 and ell[((x + s) % n, t)] == 1 and ell[((x + s + t) % n, (-s) % n)] == 1 and ell[((x + t) % n, (-t) % n)] == 1]
t4 = all(ell[(x, 2)] == 1 for x in (0, 2, 4, 6))
AB = all(ell[(x, s)] + ell[((x + s) % n, t)] == ell[(x, t)] + ell[((x + t) % n, s)] for x in range(n) for s in S for t in S)
a = {s: Fr(sum(ell[(x, s)] for x in range(n)), n) for s in S}
xi = {s: (a[s] / 4) % 1 for s in S}
chars = [j for j in range(n) if all(Fr(j * s, n) % 1 == xi[s] for s in S)]
kappa = max(min(min(Fr(j * s % n, n), 1 - Fr(j * s % n, n)) for s in S) for j in range(n))
print("colouring", [c[x] for x in range(n)], "tight 4-cycle 0-2-4-6:", t4, "tight squares:", sq, "A=B everywhere:", AB)
print("a(s):", a, "xi = a/4 equals the character x -> j x/8 for j in", chars, "kappa(S) =", kappa)
