"""max over b != 0 of |sum_{z in mu_{q+1}} exp(2 pi i Tr(bz)/p)| versus 2 sqrt(q) (the Weil-type bound for the
norm-one torus).  Uses the trace table of kappa1.py."""
import sys, cmath, math, random
sys.path.insert(0, '.')
from kappa1 import primitive_poly, polymulmod, polypow
def tables(p, f, seed=1):
    rng = random.Random(seed); n = 2 * f; q = p ** f
    m = primitive_poly(p, n, rng)
    def trace(a):
        s = [0] * n; y = a
        for _ in range(n):
            s = [(u + v) % p for u, v in zip(s, y)]; y = polypow(y, p, m, p)
        return s[0]
    T0 = []; g = [0, 1] + [0] * (n - 2); y = [1] + [0] * (n - 1)
    for k in range(n):
        T0.append(trace(y)); y = polymulmod(y, g, m, p)
    Nn = q * q - 1; T = T0 + [0] * (Nn - n)
    for k in range(n, Nn):
        T[k] = (-sum(m[i] * T[k - n + i] for i in range(n))) % p
    return q, T
for p, f in [(3, 1), (7, 1), (11, 1), (19, 1), (3, 3), (7, 3), (3, 5), (2, 1), (2, 3), (2, 5)]:
    q, T = tables(p, f)
    w = [cmath.exp(2j * math.pi * a / p) for a in range(p)]
    best = 0.0
    for j in range(q - 1):            # b = g^j, z = g^{(q-1) m}
        s = sum(w[T[j + (q - 1) * mm]] for mm in range(q + 1))
        best = max(best, abs(s))
    zero_free = None
    print(f"p={p} f={f} q={q}: max |S(b)| = {best:.4f}, 2 sqrt(q) = {2 * math.sqrt(q):.4f}, ratio {best / (2 * math.sqrt(q)):.4f};"
          f" (q+1)/p - (p-1)/p * max|S| = {(q + 1) / p - (p - 1) / p * best:.3f}", flush=True)

# The identity S(b) = -Kl(N(b)), Kl(c) = sum_{y in F_q^*} psi(Tr_{F_q/F_p}(y + c/y)), checked for every b.
def check_identity(p, f, seed=1):
    rng = random.Random(seed); n = 2 * f; q = p ** f
    m = primitive_poly(p, n, rng)
    Q2 = q * q - 1
    anti = []; y = [1] + [0] * (n - 1); g = [0, 1] + [0] * (n - 2)
    for k in range(Q2):
        anti.append(tuple(y)); y = polymulmod(y, g, m, p)
    log = {v: k for k, v in enumerate(anti)}
    q2_, T = tables(p, f, seed)
    inv2 = pow(2, -1, p)
    w = [cmath.exp(2j * math.pi * a / p) for a in range(p)]
    worst = 0.0
    for j in range(0, q - 1):                      # b = g^j; N(b) = g^{j(q+1)}
        S = sum(w[T[j + (q - 1) * mm]] for mm in range(q + 1))
        ec = (j * (q + 1)) % Q2
        Kl = 0
        for mm in range(q - 1):                    # y = g^{(q+1) mm} runs over F_q^*
            ey = ((q + 1) * mm) % Q2
            s_ = tuple((a + b) % p for a, b in zip(anti[ey], anti[(ec - ey) % Q2]))
            if any(s_):
                tr = T[log[s_]] * inv2 % p          # Tr_{F_q/F_p}(x) = Tr_{F_{q^2}/F_p}(x)/2 for x in F_q
            else:
                tr = 0
            Kl += w[tr]
        worst = max(worst, abs(S + Kl))
    print(f"identity S(b) = -Kl(N(b)) for p={p} f={f}: max |S + Kl| over all b = {worst:.2e}", flush=True)

for p, f in [(3, 1), (7, 1), (11, 1), (19, 1), (3, 3), (7, 3)]:
    check_identity(p, f)
