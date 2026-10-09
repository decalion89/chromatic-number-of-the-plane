"""Referee tests (own code) of Lemma A (tight cycles + the rounding construction of its proof) and Lemma B
(squares, winding independent of the start, averaged character, tight walks) on small examples.

Lemma A test: for random graphs on <= 7 vertices and a few named graphs, chi_c is computed as the least p/q
(p <= n, by Vince / Bondy-Hell) with a (p,q)-colouring; then for every (p,q), p <= n+2, with a colouring:
   [every (p,q)-colouring has a tight cycle]  <=>  p/q == chi_c      (Lemma A is "=>"; "<=" is Guichard's half)
and for p/q > chi_c the colouring built exactly as in the draft's proof (g = (p/r) f, c = floor(g) mod p, from an
optimal colouring f) is checked to be a (p,q)-colouring with acyclic tight digraph.
Lemma B test: Cay(Z/m, S), S = -S generating, all (p,q)-colourings with 2q <= p < 4q (c(0) = 0); checks
   (1) squares, (2) Lambda(W, x) independent of x for many closed walks, (3) a(-s) = p - a(s) and xi = a/p is a
   character of Z/m, (4) for every positive relation (small multiplicities) among T = {u : a(u) = q}, every
   ordering of the walk from every start is tight at every step; also (5) for p >= 4q the square identity can fail.
"""
import random, itertools
from fractions import Fraction as Fr
from math import floor, gcd

random.seed(12345)

# ---------------------------------------------------------------- generic colouring tools
def colourings(n, adj, p, q, fix0=True, limit=None):
    res = []
    c = [None] * n
    order = list(range(n))
    def ok(v, col):
        for u in adj[v]:
            if c[u] is not None:
                d = (col - c[u]) % p
                if not (q <= d <= p - q):
                    return False
        return True
    def rec(i):
        if limit is not None and len(res) >= limit:
            return
        if i == n:
            res.append(tuple(c))
            return
        v = order[i]
        for col in ([0] if (fix0 and i == 0) else range(p)):
            if ok(v, col):
                c[v] = col
                rec(i + 1)
                c[v] = None
    rec(0)
    return res

def tight_arcs(edges, c, p, q):
    arcs = []
    for x, y in edges:
        if (c[y] - c[x]) % p == q:
            arcs.append((x, y))
        if (c[x] - c[y]) % p == q:
            arcs.append((y, x))
    return arcs

def has_cycle(n, arcs):
    out = [[] for _ in range(n)]
    for x, y in arcs:
        out[x].append(y)
    st = [0] * n
    def dfs(v):
        st[v] = 1
        for w in out[v]:
            if st[w] == 1 or (st[w] == 0 and dfs(w)):
                return True
        st[v] = 2
        return False
    return any(st[v] == 0 and dfs(v) for v in range(n))

# ---------------------------------------------------------------- Lemma A
def lemmaA_graph(n, edges, name):
    adj = [[] for _ in range(n)]
    for x, y in edges:
        adj[x].append(y); adj[y].append(x)
    fracs = sorted({Fr(p, q) for p in range(2, n + 3) for q in range(1, p // 2 + 1)})
    feas = {}
    for fr in fracs:
        p, q = fr.numerator, fr.denominator
        feas[fr] = bool(colourings(n, adj, p, q, limit=1))
    chic = min(fr for fr in fracs if feas[fr] and fr.numerator <= n)
    # sanity: nothing smaller with p <= n+2
    assert not any(feas[fr] for fr in fracs if fr < chic)
    for fr in fracs:
        if not feas[fr]:
            continue
        for mult in (1, 2):                     # also non-reduced (p,q) = (mult*P, mult*Q) when small
            p, q = fr.numerator * mult, fr.denominator * mult
            if p > n + 2:
                continue
            cols = colourings(n, adj, p, q)
            all_tight = all(has_cycle(n, tight_arcs(edges, c, p, q)) for c in cols)
            assert all_tight == (fr == chic), (name, p, q, chic, all_tight)
            if fr > chic:
                # rounding construction of the draft's proof, from an optimal colouring
                P0, Q0 = chic.numerator, chic.denominator
                c0 = colourings(n, adj, P0, Q0, limit=1)[0]
                r = Fr(P0, Q0)
                f = [Fr(c0[v], Q0) for v in range(n)]          # in R / rZ, distance >= 1 on edges
                g = [(Fr(p) / r) * f[v] for v in range(n)]     # in R / pZ, in [0, p)
                assert all(0 <= x < p for x in g)
                c = [floor(x) % p for x in g]
                assert all(q <= (c[y] - c[x]) % p <= p - q for x, y in edges), (name, p, q)
                assert not has_cycle(n, tight_arcs(edges, c, p, q)), (name, p, q)
    return chic

def random_graph(n, pr):
    return [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < pr]

named = {
    "C5": (5, [(i, (i + 1) % 5) for i in range(5)]),
    "C7": (7, [(i, (i + 1) % 7) for i in range(7)]),
    "K4": (4, list(itertools.combinations(range(4), 2))),
    "W5 (odd wheel)": (6, [(i, (i + 1) % 5) for i in range(5)] + [(i, 5) for i in range(5)]),
    "Moser": (7, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4), (0, 5), (4, 5), (4, 6), (5, 6), (3, 6)]),
    "K7-minus-matching": (7, [e for e in itertools.combinations(range(7), 2) if e not in [(0, 1), (2, 3), (4, 5)]]),
}
countA = 0
seen_values = set()
for name, (n, E) in named.items():
    v = lemmaA_graph(n, E, name)
    seen_values.add(v)
    print(f"Lemma A on {name}: chi_c = {v}  (equivalence and rounding construction verified)")
    countA += 1
for t in range(300):
    n = random.randint(3, 7)
    E = random_graph(n, random.choice([0.3, 0.45, 0.6, 0.75]))
    if not E:
        continue
    v = lemmaA_graph(n, E, f"random{t}")
    seen_values.add(v)
    countA += 1
print(f"Lemma A: {countA} graphs passed; chi_c values seen: {sorted(seen_values)}")

# ---------------------------------------------------------------- Lemma B on Cay(Z/m, S)
def lemmaB_case(m, S, p, q, max_cols=400):
    S = sorted(set(S) | {(-s) % m for s in S})
    assert 0 not in S
    assert gcd(m, *S) == 1
    adj = [[(x + s) % m for s in S] for x in range(m)]
    adjs = [sorted(set(a)) for a in adj]
    edges = sorted({tuple(sorted((x, (x + s) % m))) for x in range(m) for s in S})
    cols = colourings(m, adjs, p, q, limit=max_cols)
    stats = {"cols": len(cols), "rel": 0, "tightcheck": 0, "square_fail": 0}
    for c in cols:
        def ell(x, s):
            d = (c[(x + s) % m] - c[x]) % p
            assert q <= d <= p - q
            return d
        # (1) squares
        for x in range(m):
            for s in S:
                for t in S:
                    lhs = ell(x, s) + ell((x + s) % m, t)
                    rhs = ell(x, t) + ell((x + t) % m, s)
                    if p < 4 * q:
                        assert lhs == rhs
                    elif lhs != rhs:
                        stats["square_fail"] += 1
        if p >= 4 * q:
            continue
        # (2) winding constant for random closed walks
        for _ in range(30):
            L = random.randint(2, 8)
            steps = [random.choice(S) for _ in range(L - 1)]
            last = (-sum(steps)) % m
            if last == 0 or last not in S:
                continue
            steps.append(last)
            vals = set()
            for x in range(m):
                y, tot = x, 0
                for s in steps:
                    tot += ell(y, s); y = (y + s) % m
                assert y == x and tot % p == 0
                vals.add(tot)
            assert len(vals) == 1
        # (3) averaging (the invariant mean on Z/m is the plain average)
        a = {s: Fr(sum(ell(x, s) for x in range(m)), m) for s in S}
        for s in S:
            assert q <= a[s] <= p - q and a[(-s) % m] == p - a[s]
        ok = [j for j in range(m) if all((a[s] / p - Fr(j * s, m)).denominator == 1 for s in S)]
        assert len(ok) == 1, "xi is not a character"
        # (4) tight walks
        T = [u for u in S if a[u] == q]
        for r in (1, 2, 3):
            for us in itertools.combinations(T, r):
                for ns in itertools.product(range(1, 5), repeat=r):
                    if sum(nn * u for nn, u in zip(ns, us)) % m != 0:
                        continue
                    stats["rel"] += 1
                    walk = [u for nn, u in zip(ns, us) for _ in range(nn)]
                    perms = set(itertools.permutations(walk)) if len(walk) <= 7 else {tuple(random.sample(walk, len(walk))) for _ in range(40)}
                    for w in perms:
                        for x in range(m):
                            y = x
                            for s in w:
                                assert ell(y, s) == q, ("tight walk not tight", m, S, p, q, c, w, x)
                                y = (y + s) % m
                            stats["tightcheck"] += 1
                    # and the tight digraph of the colouring has a directed cycle
                    assert has_cycle(m, tight_arcs(edges, c, p, q))
    return stats

casesB = [
    (7, [1, 2], 3, 1), (9, [1, 3], 3, 1), (12, [1, 4], 3, 1), (8, [1, 3], 3, 1), (10, [1, 3], 3, 1),
    (13, [1, 5], 3, 1), (11, [1, 3], 5, 2), (12, [1, 5], 5, 2), (14, [1, 3], 7, 2), (13, [2, 5], 7, 2),
    (15, [1, 4], 7, 2), (16, [1, 6], 7, 2), (17, [1, 4], 7, 3), (12, [1, 5], 8, 3), (18, [1, 5], 7, 2),
    (21, [1, 6], 7, 2), (20, [3, 7], 7, 2),
]
tot_rel = tot_tc = tot_cols = 0
for (m, S, p, q) in casesB:
    st = lemmaB_case(m, S, p, q)
    tot_rel += st["rel"]; tot_tc += st["tightcheck"]; tot_cols += st["cols"]
    print(f"Lemma B on Cay(Z/{m}, +-{S}), (p,q)=({p},{q}): {st['cols']} colourings; "
          f"{st['rel']} positive tight relations, {st['tightcheck']} (walk, start) pairs all tight")
print(f"Lemma B: {tot_cols} colourings, {tot_rel} relations, {tot_tc} tight walk checks passed")
# (5) p >= 4q: the square identity can fail
fails = 0
for (m, S, p, q) in [(10, [1, 3], 8, 2), (11, [1, 4], 9, 2), (13, [1, 5], 4, 1), (9, [1, 2], 4, 1)]:
    st = lemmaB_case(m, S, p, q, max_cols=200)
    fails += st["square_fail"]
    print(f"p >= 4q control: Cay(Z/{m}, +-{S}), ({p},{q}): {st['cols']} colourings, square-identity failures {st['square_fail']}")
print("controls: square identity fails for some colouring with p >= 4q:", fails > 0)
