"""hoffman_padic.py PMAX: for every prime p = 3 (mod 4) below PMAX, the eigenvalues of the unit-distance graph of
F_p^2 (the p + 1 unit vectors), one per norm class of characters, and Hoffman's ratio -mu/(1 - mu), mu = the
smallest eigenvalue over the degree p + 1. A ratio below 1/4 means no proper 4-colouring; since the eigenvalues
of the higher levels (Z/p^k)^2 are at most 2/(p + 1) of the degree (spectrum.py), the same ratio bounds every
level. Also checks |lambda| <= 2 sqrt(p)."""
import sys
import numpy as np
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 104


def isprime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))


for p in range(3, PMAX, 4):
    if not isprime(p):
        continue
    C = np.array([(x, y) for x in range(p) for y in range(p) if (x * x + y * y) % p == 1])
    assert len(C) == p + 1
    reps = {}
    for x in range(p):
        for y in range(p):
            N = (x * x + y * y) % p
            if N and N not in reps:
                reps[N] = (x, y)
    X = np.array(list(reps.values()))
    lam = np.cos(2 * np.pi * ((X[:, :1] * C[None, :, 0] + X[:, 1:] * C[None, :, 1]) % p) / p).sum(axis=1)
    assert np.abs(lam).max() <= 2 * np.sqrt(p) + 1e-9
    mu = lam.min() / (p + 1)
    ratio = -mu / (1 - mu)
    print(f"p = {p}: smallest eigenvalue {lam.min():.10f}, Hoffman ratio {ratio:.7f}"
          f" ({'below' if ratio < 0.25 else 'not below'} 1/4)", flush=True)
