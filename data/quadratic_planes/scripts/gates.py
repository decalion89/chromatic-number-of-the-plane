"""Valuations of the unit vectors z = ((a + b r) + i (c + e r)) / D, r = sqrt(d), at the places of K(i) above p,
for the 'gates' (places where every direction is integral give a homomorphism to a finite graph).
Split case (p splits in K(i)/K): v_P(z) = -v_cP(z); we report the multiset of v_P(z).
We compute v_P through the norm to a quadratic subfield in which the place below P is inert or ramified."""
import sys
from math import gcd
from units_fast import units_fast


def vp(n, p):
    if n == 0:
        return 99
    k = 0
    while n % p == 0:
        n //= p; k += 1
    return k


def padic_sqrt(a, p, k):
    """root t of t^2 = a mod p^k, a a nonzero square mod p (p odd)"""
    m = p
    t = next(x for x in range(p) if (x * x - a) % p == 0)
    while m < p ** k:
        m2 = m * p
        # Newton step
        t = (t - (t * t - a) * pow(2 * t, -1, m2)) % m2
        m = m2
    return t


def val_split_gauss(a, b, c, e, D, d, p):
    """p = 1 mod 4: Gaussian prime pi | p; place above pi in K(i). Uses N_{K(i)/Q(i)}(z) = ((a+ic)^2 - d (b+ie)^2)/D^2.
    Valid when p ramifies in K (p | d): then v_P(z) = v_pi(N z)."""
    k = 30
    iota = padic_sqrt(-1, p, k)             # image of i in Z_p for the prime pi = (p, i - iota)
    M = p ** k
    x = (a + iota * c) % M
    y = (b + iota * e) % M
    N = (x * x - d * y * y) % M
    return vp(N, p) - 2 * vp(D, p) if N else 99


def val_via_minus_d(a, b, c, e, D, d, p, sign=1):
    """-d a nonzero square mod p: the prime q | p of Q(sqrt(-d)) with sqrt(-d) -> sign*t; P | q in K(i) inert over q.
    N_{K(i)/Q(sqrt-d)}(z) = (X + Y sqrt(-d))/D^2 with X = a^2 + c^2 - d (b^2 + e^2), Y = 2 (a e - b c).
    v_P(z) = (v_q(X + Y t) - 2 v_p(D)) / 2."""
    k = 30
    t = sign * padic_sqrt(-d % p ** k, p, k)
    M = p ** k
    X = a * a + c * c - d * (b * b + e * e)
    Y = 2 * (a * e - b * c)
    N = (X + Y * t) % M
    v = vp(N, p) if N else 99
    return (v - 2 * vp(D, p)) / 2


if __name__ == "__main__":
  d = int(sys.argv[1]); Ds = [int(x) for x in sys.argv[2].split(",")]
  for D in Ds:
      U = units_fast(d, D)
      out = [f"d={d} D={D}: {len(U)} units"]
      for p in (3, 5, 7, 11, 13, 17, 19, 23, 29):
          if d % p == 0 and p % 4 == 1:
              vs = [val_split_gauss(*u, D, d, p) for u in U]
              out.append(f"p={p} (ram in K, split in K(i)): v_P != 0 for {sum(1 for v in vs if v != 0)}")
          elif d % p and pow(-d % p, (p - 1) // 2, p) == 1 and pow(d % p, (p - 1) // 2, p) != 1:
              # -d square, d non-square: p inert in K, splits in K(i) (residue F_{p^2}) -> gate H_{p^2}
              vs = [val_via_minus_d(*u, D, d, p) for u in U]
              out.append(f"p={p} (inert in K, split in K(i)): v_P != 0 for {sum(1 for v in vs if v != 0)}")
      print("; ".join(out))


def val_split_complete(a, b, c, e, D, d, p, s_sign=1, i_sign=1):
    """p = 1 mod 4 and d a nonzero square mod p: K(i) embeds in Q_p; place given by sqrt(d) -> s_sign*s, i -> i_sign*iota."""
    k = 30
    M = p ** k
    s = s_sign * padic_sqrt(d % M, p, k)
    iota = i_sign * padic_sqrt(-1 % M, p, k)
    N = ((a + b * s) + iota * (c + e * s)) % M
    return (vp(N, p) if N else 99) - vp(D, p)


if __name__ == "__main__" and len(sys.argv) > 3:
    for D in Ds:
        U = units_fast(d, D)
        for p in (5, 13, 17, 29):
            if p % 4 == 1 and d % p and pow(d % p, (p - 1) // 2, p) == 1:
                cnt = sum(1 for u in U if val_split_complete(*u, D, d, p) != 0)
                print(f"d={d} D={D}: p={p} split completely in K(i) (H_{p}): v_P != 0 for {cnt} of {len(U)}")
