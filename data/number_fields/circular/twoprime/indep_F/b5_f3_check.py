"""Referee: exhaustive check of Lemma B5' for p = 7, f = 3 (and f = 1), by direct computation.
F_{7^6} = F_49[t]/(t^3 - a) with a a non-cube of F_49; Tr_{F_{7^6}/F_49}(c0 + c1 t + c2 t^2) = 3 c0.
Every F_49-linear map F_{7^6} -> F_49 is x -> Tr(b x).  Claim: no b maps mu_{344} into A'_7.
(Also: for f = 1, x -> e x with e in A'_7 maps mu_8 into A'_7.)"""
P = 7


def m49(u, v):
    return ((u[0] * v[0] - u[1] * v[1]) % P, (u[0] * v[1] + u[1] * v[0]) % P)


def a49(u, v):
    return ((u[0] + v[0]) % P, (u[1] + v[1]) % P)


def conj(u):
    return (u[0], (-u[1]) % P)


F49 = [(x, y) for x in range(P) for y in range(P)]
nz = [u for u in F49 if u != (0, 0)]
mu8 = [z for z in F49 if m49(z, conj(z)) == (1, 0)]
Ap = set(e for e in F49 if all(m49(conj(e), z)[0] in (2, 3, 4, 5) for z in mu8))
assert len(Ap) == 8


def pw49(u, n):
    r = (1, 0)
    for _ in range(n):
        r = m49(r, u)
    return r


cubes = set(pw49(u, 3) for u in nz)
a = next(u for u in nz if u not in cubes)
print("F_49 non-cube a =", a, "; number of cubes:", len(cubes))
ZERO, ONE = (0, 0), (1, 0)


def mul(X, Y):
    # (x0 + x1 t + x2 t^2)(y0 + y1 t + y2 t^2), t^3 = a
    c = [ZERO] * 5
    for i in range(3):
        for j in range(3):
            c[i + j] = a49(c[i + j], m49(X[i], Y[j]))
    return (a49(c[0], m49(a, c[3])), a49(c[1], m49(a, c[4])), c[2])


def pw(X, n):
    R, B = (ONE, ZERO, ZERO), X
    while n:
        if n & 1:
            R = mul(R, B)
        B = mul(B, B)
        n >>= 1
    return R


Q2 = 7 ** 6
q = 7 ** 3
order = Q2 - 1
n_ = order
pf = []
d = 2
while d * d <= n_:
    if n_ % d == 0:
        pf.append(d)
        while n_ % d == 0:
            n_ //= d
    d += 1
if n_ > 1:
    pf.append(n_)
print("7^6 - 1 =", order, "prime factors", pf)
elts = [(x, y, z) for x in F49 for y in F49 for z in F49]
g = None
for X in elts:
    if X == (ZERO, ZERO, ZERO):
        continue
    if all(pw(X, order // p_) != (ONE, ZERO, ZERO) for p_ in pf):
        g = X
        break
print("generator found:", g is not None)
zeta = pw(g, q - 1)            # generator of mu_{q+1}
mu = [(ONE, ZERO, ZERO)]
for _ in range(q + 1):
    mu.append(mul(mu[-1], zeta))
assert mu[-1] == (ONE, ZERO, ZERO) and len(set(mu[:-1])) == q + 1
mu = mu[:-1]
# norm-one check: z^(q+1) = 1
assert all(pw(z, q + 1) == (ONE, ZERO, ZERO) for z in mu[:5])
three = (3, 0)
bad = 0
b = (ONE, ZERO, ZERO)
gpow = (ONE, ZERO, ZERO)
nb = 0
for beta in range(q - 1):          # b = g^beta runs over representatives of F^x / mu_{q+1}
    ok = True
    for z in mu:
        val = m49(three, mul(gpow, z)[0])
        if val not in Ap:
            ok = False
            break
    if ok:
        bad += 1
    nb += 1
    gpow = mul(gpow, g)
print(f"f = 3: checked {nb} classes of b (b != 0, modulo mu_{q+1}); maps with lambda(mu_{q+1}) in A'_7: {bad}")
# f = 1: lambda(x) = e x
print("f = 1: e*mu_8 subset of A'_7 for e in A'_7:", all(set(m49(e, z) for z in mu8) <= Ap for e in Ap))
