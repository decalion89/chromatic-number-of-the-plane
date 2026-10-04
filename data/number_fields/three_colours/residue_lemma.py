"""Local lemma at 3 for the three-colour local-global theorem.

For q = 3^f (f odd), F_{q^2} contains F_9 and mu_{q+1}.  Claim: for f >= 3 there is no b in F_{q^2} with
Tr_{F_{q^2}/F_9}(b*z) in A = {x in F_9 : x^4 = -1} (= (1+i)mu_4) for every z in mu_{q+1};
for f = 1 such b exist (b in A).  Exhaustive over b, exact arithmetic in F_{3^{2f}} = F_3[t]/(irreducible)."""
import itertools, sys

def poly_irreducible(n, p=3):
    # find a primitive polynomial of degree n over F_p by brute force (monic), test primitivity via order
    for coeffs in itertools.product(range(p), repeat=n):
        if coeffs[0] == 0: continue
        f = list(coeffs) + [1]          # f = c0 + c1 t + ... + t^n
        # build field via log table attempt: check t generates a cyclic group of order p^n - 1
        q = p ** n
        # multiplication by t on vectors
        def mul_t(v):
            carry = v[-1]
            w = [0] + v[:-1]
            return [(w[k] - carry * f[k]) % p for k in range(n)]
        v = [1] + [0] * (n - 1)
        seen = 0
        order = None
        for e in range(1, q):
            v = mul_t(v)
            if v == [1] + [0] * (n - 1):
                order = e; break
        if order == q - 1:
            return f
    raise ValueError

def build(n, p=3):
    f = poly_irreducible(n, p)
    q = p ** n
    def enc(v): return sum(c * p ** k for k, c in enumerate(v))
    exp = [0] * (q - 1); log = {}
    v = [1] + [0] * (n - 1)
    for e in range(q - 1):
        exp[e] = enc(v); log[enc(v)] = e
        carry = v[-1]; w = [0] + v[:-1]
        v = [(w[k] - carry * f[k]) % p for k in range(n)]
    def dec(x): return [(x // p ** k) % p for k in range(n)]
    def add(x, y): return enc([(a + b) % p for a, b in zip(dec(x), dec(y))])
    return q, exp, log, add

def check(f):
    n = 2 * f
    Q, exp, log, add = build(n)
    q = 3 ** f
    g = Q - 1                      # multiplicative group order
    # F_9 inside: elements x with x^9 = x: exponents multiple of (Q-1)/8
    s9 = (Q - 1) // 8
    F9 = [0] + [exp[(k * s9) % g] for k in range(8)]
    # A = elements of F_9 with x^4 = -1: exponents e with 4e = (Q-1)/2 mod Q-1 -> e = s9*(odd)
    A = set(exp[(k * s9) % g] for k in (1, 3, 5, 7))
    mu = [exp[(k * (q - 1)) % g] for k in range(q + 1)]   # mu_{q+1}: exponents multiples of q-1
    def mul(x, y):
        if x == 0 or y == 0: return 0
        return exp[(log[x] + log[y]) % g]
    def frob9(x):   # x -> x^9
        return 0 if x == 0 else exp[(9 * log[x]) % g]
    def tr(x):      # Tr_{F_{Q}/F_9}(x) = sum_{j<f} x^{9^j}
        s, y = 0, x
        for _ in range(f):
            s = add(s, y); y = frob9(y)
        return s
    good = []
    # b up to the action of mu_{q+1} * F_9^x (the condition is invariant): just scan all b
    for b in range(1, Q):
        ok = True
        for z in mu:
            if tr(mul(b, z)) not in A:
                ok = False; break
        if ok:
            good.append(b)
            if len(good) > 3: break
    return good

for f in [1, 3, 5]:
    good = check(f)
    print(f"f = {f}: q = {3**f}: b with Tr(b z) in A for all z in mu_(q+1): {'none' if not good else str(len(good))+'+ found'}", flush=True)
