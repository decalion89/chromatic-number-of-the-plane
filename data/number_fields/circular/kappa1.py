"""kappa_1(p, f): the best level-1 (residue-field) circular colouring at a place above p with residue field F_q,
q = p^f, i not in F_v.  kappa_1 = max over F_p-linear lam: F_{q^2} -> F_p of min over mu_{q+1} of ||lam(z)/p||.
Every lam is y -> Tr(b y); b and b*zeta (zeta in mu_{q+1}) give the same value set, so b runs over g^j, 0 <= j < q-1.
Tr(g^k) is a linear recurrence with the minimal polynomial of g."""
import sys, random
from fractions import Fraction

def polymulmod(a, b, m, p):
    n = len(m) - 1
    res = [0] * (2 * n)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                res[i + j] = (res[i + j] + x * y) % p
    for k in range(2 * n - 1, n - 1, -1):
        c = res[k]
        if c:
            for i in range(n + 1):
                res[k - n + i] = (res[k - n + i] - c * m[i]) % p
    return res[:n]

def polypow(a, e, m, p):
    n = len(m) - 1
    r = [1] + [0] * (n - 1)
    while e:
        if e & 1:
            r = polymulmod(r, a, m, p)
        a = polymulmod(a, a, m, p)
        e >>= 1
    return r

def factor(n):
    fs, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            fs.add(d); n //= d
        d += 1
    if n > 1:
        fs.add(n)
    return fs

def primitive_poly(p, n, rng):
    Q = p ** n - 1
    fs = factor(Q)
    one = [1] + [0] * (n - 1)
    x = [0, 1] + [0] * (n - 2)
    while True:
        m = [rng.randrange(p) for _ in range(n)] + [1]
        if m[0] == 0:
            continue
        if polypow(x, Q, m, p) != one:
            continue
        if all(polypow(x, Q // r, m, p) != one for r in fs):
            return m

def kappa1(p, f, seed=1):
    rng = random.Random(seed)
    n = 2 * f
    q = p ** f
    m = primitive_poly(p, n, rng)
    # traces of g^k for k < n by explicit Frobenius sums
    def trace(a):
        s = [0] * n
        y = a
        for _ in range(n):
            s = [(u + v) % p for u, v in zip(s, y)]
            y = polypow(y, p, m, p)
        assert all(c == 0 for c in s[1:])
        return s[0]
    T0 = []
    g = [0, 1] + [0] * (n - 2)
    y = [1] + [0] * (n - 1)
    for k in range(n):
        T0.append(trace(y))
        y = polymulmod(y, g, m, p)
    N = q * q - 1
    T = T0 + [0] * (N - n)
    for k in range(n, N):
        s = 0
        for i in range(n):
            s -= m[i] * T[k - n + i]
        T[k] = s % p
    best, bestj = -1, None
    for j in range(q - 1):
        mn = p
        for mm in range(q + 1):
            t = T[j + (q - 1) * mm]
            d = min(t, p - t)
            if d < mn:
                mn = d
                if mn <= best:
                    break
        if mn > best:
            best, bestj = mn, j
    return Fraction(best, p), bestj

if __name__ == "__main__":
    for p, f in [(int(a), int(b)) for a, b in (s.split(",") for s in sys.argv[1:])]:
        k, j = kappa1(p, f)
        print(f"p={p} f={f} q={p**f}: kappa_1 = {k} = {float(k):.4f}  (> 1/4: {k > Fraction(1, 4)})  b = g^{j}", flush=True)
