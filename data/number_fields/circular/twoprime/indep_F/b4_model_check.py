"""Referee: finite-model check of Lemma B4' (level one at 7) for residue degree f = 1 and ramification e = 1, 2, 3.
O/7O = F_49[pi]/(pi^e), sigma(sum a_k pi^k) = sum a_k^7 pi^k (pi in F_u is fixed, sigma = Frobenius on Z_49).
The image of T(F_u) in O/7O is the set of z with z sigma(z) = 1 (norm-one elements; the norm is onto on higher
unit groups of an unramified extension, so every such z lifts).  An F_49-linear lambda: O/7O -> F_49 is
lambda(sum a_k pi^k) = sum c_k a_k.  Lemma B4' predicts: lambda(T) in A'_7 forces c_k = 0 for k >= 1 (lambda kills
the maximal ideal), and then c_0 mu_8 in A'_7."""
from itertools import product
P = 7


def m49(u, v):
    return ((u[0] * v[0] - u[1] * v[1]) % P, (u[0] * v[1] + u[1] * v[0]) % P)


def a49(u, v):
    return ((u[0] + v[0]) % P, (u[1] + v[1]) % P)


def conj(u):
    return (u[0], (-u[1]) % P)


F49 = [(x, y) for x in range(P) for y in range(P)]
mu8 = [z for z in F49 if m49(z, conj(z)) == (1, 0)]
Ap = set(e for e in F49 if all(m49(conj(e), z)[0] in (2, 3, 4, 5) for z in mu8))


def frob(u):  # u^7 = conjugate in F_49
    return conj(u)


for e in (1, 2, 3):
    def mul(X, Y):
        c = [(0, 0)] * e
        for i in range(e):
            for j in range(e - i):
                c[i + j] = a49(c[i + j], m49(X[i], Y[j]))
        return tuple(c)
    one = tuple([(1, 0)] + [(0, 0)] * (e - 1))
    T = [z for z in product(F49, repeat=e) if mul(z, tuple(frob(a) for a in z)) == one]
    expected = 8 * 7 ** (e - 1)
    good = []
    for c in product(F49, repeat=e):
        ok = True
        for z in T:
            v = (0, 0)
            for k in range(e):
                v = a49(v, m49(c[k], z[k]))
            if v not in Ap:
                ok = False
                break
        if ok:
            good.append(c)
    kill = all(all(ck == (0, 0) for ck in c[1:]) for c in good)
    c0s = sorted(set(c[0] for c in good))
    print(f"e = {e}: |T mod 7O| = {len(T)} (expected 8*7^(e-1) = {expected}); maps with lambda(T) in A'_7: "
          f"{len(good)}; all kill the maximal ideal: {kill}; their c_0 = A'_7: {c0s == sorted(Ap)}")
