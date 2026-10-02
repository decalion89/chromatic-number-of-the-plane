"""padic_reach.py: how far the certified graphs reach among the p-adic planes (note §5).

chi(Q_p^2) >= 4 as soon as one of the fields Q(sqrt d) of data/quadratic_planes/ embeds in Q_p, that is when d is a
nonzero square mod p (p odd; Hensel): each q{d}.json, d = 47 included, is a unit-distance graph with no 3-colouring.
For p = 1 (mod 4) the plane needs infinitely many colours anyway (Davies; note §5), so the question is p = 3 (mod 4).
No finite set of fields reaches every such p (Chebotarev; note §5). padic_reach.c, a segmented sieve with Jacobi
symbols, finds the first prime p = 3 (mod 4) after 3 that none of the 27 fields reaches, P0 = 2 129 503 819 (21 s).
So every prime p = 3 (mod 4) with 7 <= p < P0 has some d as a nonzero square, and chi(Q_p^2) >= 4 for all of them.

This script checks, independently of the C program: that its list of d is the list of graphs in the data; that P0
is prime (deterministic Miller-Rabin), is 3 mod 4, and has no d as a nonzero square (Euler's criterion); and, by
Euler's criterion, that every prime p = 3 (mod 4) with 7 <= p < N (default 10^7) has one.
usage: python3 padic_reach.py [N]"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
P0 = 2129503819


def fields():
    return sorted(int(f[1:-5]) for f in os.listdir(DATA) if re.fullmatch(r"q\d+\.json", f))


def c_list():
    src = open(os.path.join(HERE, "padic_reach.c"), encoding="utf-8").read()
    body = re.search(r"static const int D\[\] = \{([^}]*)\}", src).group(1)
    return sorted(int(x) for x in body.replace("\n", " ").split(",") if x.strip())


def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)   # deterministic for n < 3.3e24
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def reached(p, ds):
    """the d that are nonzero squares mod the odd prime p"""
    return [d for d in ds if d % p and pow(d, (p - 1) // 2, p) == 1]


def check(N=10 ** 7):
    ds = fields()
    assert len(ds) == 27 and 47 in ds and c_list() == ds, "padic_reach.c and the data disagree"
    assert is_prime(P0) and P0 % 4 == 3 and reached(P0, ds) == []
    sieve = bytearray([1]) * N
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, N, i)))
    missed = [p for p in range(3, N, 4) if sieve[p] and not reached(p, ds)]
    assert missed == [3], missed
    return ds, sum(1 for p in range(3, N, 4) if sieve[p])


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 10 ** 7
    ds, count = check(N)
    print(f"{len(ds)} fields; P0 = {P0} is prime, = 3 mod 4, and no d is a nonzero square mod P0")
    print(f"every one of the {count} primes p = 3 (mod 4) below {N} except p = 3 has some d as a nonzero square")
    print(f"so chi(Q_p^2) >= 4 for every prime p = 3 (mod 4) with 7 <= p < {N} (padic_reach.c: up to P0)")
