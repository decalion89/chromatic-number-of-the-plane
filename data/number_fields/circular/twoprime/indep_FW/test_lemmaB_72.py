"""Extra Lemma B(4) tests for (p,q) = (7,2) and (3,1) with larger multiplicities, on Cay(Z/m, S) and on
Cay(Z/m1 x Z/m2, S); all colourings (c(0)=0) up to a limit.  Same checks as test_lemmas_AB.py, part (4)."""
import itertools, random
from fractions import Fraction as Fr
random.seed(7)

def run(G, add, S, p, q, limit=3000, maxmult=14):
    idx = {g: k for k, g in enumerate(G)}
    n = len(G)
    S = list(S)
    nb = [[idx[add(g, s)] for s in S] for g in G]
    order = list(range(n))
    cols = []
    c = [None] * n
    def ok(v, col):
        for k, s in enumerate(S):
            w = nb[v][k]
            if c[w] is not None:
                d = (c[w] - col) % p
                if not (q <= d <= p - q):
                    return False
        return True
    def rec(i):
        if len(cols) >= limit:
            return
        if i == n:
            cols.append(tuple(c)); return
        v = order[i]
        for col in ([0] if i == 0 else range(p)):
            if ok(v, col):
                c[v] = col; rec(i + 1); c[v] = None
    rec(0)
    nrel = nwalk = 0
    for col in cols:
        ell = lambda v, k: (col[nb[v][k]] - col[v]) % p
        a = [Fr(sum(ell(v, k) for v in range(n)), n) for k in range(len(S))]
        T = [k for k in range(len(S)) if a[k] == q]
        for r in (1, 2, 3):
            for ks in itertools.combinations(T, r):
                rng = range(1, maxmult + 1) if r <= 2 else range(1, 8)
                for ns in itertools.product(rng, repeat=r):
                    tot = G[0]
                    for nn, k in zip(ns, ks):
                        for _ in range(nn):
                            tot = add(tot, S[k])
                    if tot != G[0]:
                        continue
                    nrel += 1
                    walk = [k for nn, k in zip(ns, ks) for _ in range(nn)]
                    for _ in range(25):
                        random.shuffle(walk)
                        for v0 in range(n):
                            v = v0
                            for k in walk:
                                assert ell(v, k) == q, "Lemma B(4) violated"
                                v = nb[v][k]
                            assert v == v0
                            nwalk += 1
    return len(cols), nrel, nwalk

cyc = lambda m: (list(range(m)), lambda a, b: (a + b) % m)
def sym(S, neg):
    out = []
    for s in S:
        for t in (s, neg(s)):
            if t not in out:
                out.append(t)
    return out

tests = []
for m, S, p, q in [(14, [4, 5, 6], 7, 2), (21, [6, 7, 8, 9], 7, 2), (7, [2, 3], 7, 2), (14, [4, 6], 7, 2),
                   (28, [8, 9, 11, 12], 7, 2), (35, [10, 11, 13, 15], 7, 2), (9, [3, 4], 3, 1), (12, [4, 5], 3, 1),
                   (15, [5, 7], 3, 1), (21, [7, 8, 9], 7, 2)]:
    G, add = cyc(m)
    Ss = sym(S, lambda s: (-s) % m)
    tests.append((f"Z/{m} S=+-{S}", G, add, Ss, p, q))
# a two-dimensional group: Z/7 x Z/7 with a 7-adic-like set
G2 = [(a, b) for a in range(7) for b in range(7)]
add2 = lambda x, y: ((x[0] + y[0]) % 7, (x[1] + y[1]) % 7)
S2 = sym([(1, 0), (0, 1), (1, 1), (1, 3)], lambda s: ((-s[0]) % 7, (-s[1]) % 7))
tests.append(("Z/7xZ/7 S=+-{(1,0),(0,1),(1,1),(1,3)}", G2, add2, S2, 7, 2))
G3 = [(a, b) for a in range(3) for b in range(6)]
add3 = lambda x, y: ((x[0] + y[0]) % 3, (x[1] + y[1]) % 6)
S3 = sym([(1, 0), (0, 2), (1, 2), (1, 1)], lambda s: ((-s[0]) % 3, (-s[1]) % 6))
tests.append(("Z/3xZ/6 S=+-{(1,0),(0,2),(1,2),(1,1)}", G3, add3, S3, 3, 1))
TOT = [0, 0, 0]
for name, G, add, S, p, q in tests:
    r = run(G, add, S, p, q)
    TOT = [x + y for x, y in zip(TOT, r)]
    print(f"{name}, (p,q)=({p},{q}): {r[0]} colourings, {r[1]} positive tight relations, {r[2]} (shuffled walk, start) checks: all tight")
print("TOTAL", TOT, "-> Lemma B(4) never violated")
