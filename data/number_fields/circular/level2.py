"""Check the level lemma on finite models: at a place above an odd prime p with L_w/F_v unramified, a character of
O_w/p^2 that is not trivial on pO_w has min over T of ||xi|| <= 1/(2p); so the best level-2 colouring is the level-1
one.  Models: Q_7 (O = Z_7[i]) and Q_27 (F_v unramified of degree 3 over Q_3, O_w = W(F_729))."""
import numpy as np, itertools, sys
from fractions import Fraction

def q7():
    M = 49
    pts = [(a, b) for a in range(M) for b in range(M) if (a * a + b * b - 1) % M == 0]
    A = np.array(pts, dtype=np.int64)                      # 56 x 2
    c = np.array([(c1, c2) for c1 in range(M) for c2 in range(M)], dtype=np.int64)
    vals = (c @ A.T) % M
    d = np.minimum(vals, M - vals).min(axis=1)
    lvl1 = (c % 7 == 0).all(axis=1)
    print(f"Q_7: |T mod 49| = {len(pts)}; best overall {Fraction(int(d.max()), M)}; "
          f"best with c not = 0 mod 7: {Fraction(int(d[~lvl1].max()), M)} (bound 1/14 = {1/14:.4f})")

def q27():
    p, n, M = 3, 6, 9
    # primitive polynomial of degree 6 over F_3, lifted to Z/9; Frobenius phi(x) = root of g congruent to x^3
    sys.path.insert(0, '.')
    from kappa1 import primitive_poly, polymulmod, polypow
    import random
    m = primitive_poly(3, 6, random.Random(5))
    g = [c % M for c in m]                                  # monic lift
    def mul(a, b):
        return polymulmod(a, b, g, M)
    def ev(poly_coeffs_of_g, y):                            # evaluate g at y (element of O/9)
        r = [0] * n
        pw = [1] + [0] * (n - 1)
        for c in poly_coeffs_of_g:
            r = [(u + c * v) % M for u, v in zip(r, pw)]
            pw = mul(pw, y)
        return r
    x = [0, 1, 0, 0, 0, 0]
    y = polypow(x, 3, g, M)
    # Newton step: y <- y - g(y)/g'(y); g'(y) is a unit; one step gives a root mod 9
    gp = [(i * g[i]) % M for i in range(1, n + 1)]
    gy, gpy = ev(g, y), ev(gp, y)
    # invert gpy in O/9 by brute force over Teichmueller-free method: gpy^(|units|-1)
    units_order = (3 ** 6 - 1) * 3 ** 6
    inv = polypow(gpy, units_order - 1, g, M)
    assert mul(inv, gpy) == [1, 0, 0, 0, 0, 0]
    corr = mul(gy, inv)
    phix = [(u - v) % M for u, v in zip(y, corr)]
    assert ev(g, phix) == [0] * n
    def phi(a):
        r = [0] * n
        pw = [1] + [0] * (n - 1)
        for c in a:
            r = [(u + c * v) % M for u, v in zip(r, pw)]
            pw = mul(pw, phix)
        return r
    sig = lambda a: phi(phi(phi(a)))
    # all elements of O/9 as integer vectors
    allel = np.array(list(itertools.product(range(M), repeat=n)), dtype=np.int64)
    # multiplication matrices: element a acts on basis; compute z*sigma(z) vectorised via structure constants
    basis = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    mult = np.zeros((n, n, n), dtype=np.int64)               # x^i * x^j = sum_k mult[i,j,k] x^k
    for i in range(n):
        for j in range(n):
            mult[i, j] = mul(basis[i], basis[j])
    S = np.array([sig(basis[i]) for i in range(n)], dtype=np.int64)   # sigma as a matrix on coefficient vectors
    sz = (allel @ S) % M
    prod = np.einsum('ai,aj,ijk->ak', allel, sz, mult) % M
    one = np.zeros(n, dtype=np.int64); one[0] = 1
    T = allel[(prod == one).all(axis=1)]
    print(f"Q_27: |T mod 9| = {len(T)} (expected 28*27 = 756)")
    # trace form: Tr(x^i x^j) = Tr of basis products; Tr(y) = sum_k phi^k(y)
    def tr(a):
        s = [0] * n
        y = a
        for _ in range(n):
            s = [(u + v) % M for u, v in zip(s, y)]
            y = phi(y)
        assert all(c == 0 for c in s[1:]), s
        return s[0]
    trb = [tr(basis[i]) for i in range(n)]
    TM = np.einsum('ijk,k->ij', mult, np.array(trb)) % M      # Tr(x^i x^j)
    W = (T @ TM) % M                                          # row z -> functional b -> Tr(bz)
    best_all, best_gen = 0, 0
    for start in range(0, len(allel), 59049):
        B = allel[start:start + 59049]
        vals = (B @ W.T) % M
        d = np.minimum(vals, M - vals).min(axis=1)
        gen = ~((B % 3) == 0).all(axis=1)
        best_all = max(best_all, int(d.max()))
        if gen.any():
            best_gen = max(best_gen, int(d[gen].max()))
    print(f"Q_27: best overall at level 2 = {Fraction(best_all, M)}; best genuine level 2 = {Fraction(best_gen, M)} "
          f"(bound 1/6)")

q7()
q27()
