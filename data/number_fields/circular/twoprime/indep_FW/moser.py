"""Referee check (own code): the Moser spindle as a unit-distance graph over Q(sqrt3, sqrt11), and chi_c = 7/2.

Field elements: dicts over the basis 1, r3, r11, r33 (r3^2 = 3, r11^2 = 11, r3 r11 = r33) with Fraction coefficients.
Embedding: A = 0, far tips D = (sqrt11/2, 1/2), G = (sqrt11/2, -1/2) (|D| = |G| = sqrt3, |D - G| = 1);
B, C = D/2 +- n/2 with n the unit normal of D; E, F likewise for G.
Checks:
 (1) all 21 pairwise squared distances computed exactly; the pairs at distance 1 are exactly the 11 spindle edges;
 (2) chi = 4 (no proper 3-colouring by exhaustion), alpha = 2;
 (3) for every p <= 14 and q with 2 <= p/q <= 4: whether a (p,q)-colouring exists (exhaustive backtracking);
     least such p/q with p <= 14 is 7/2; this is chi_c because chi_c = p/q with p <= |V| = 7 (Bondy-Hell / Vince),
     and independently chi_c >= chi_f >= n/alpha = 7/2;
 (4) the draft's colouring (0,2,4,6,3,5,1 in the order A,B,C,D,E,F,G) is a (7,2)-colouring;
 (5) every (7,2)-colouring has a tight cycle (all colourings enumerated), as Lemma A requires for chi_c = 7/2.
"""
from fractions import Fraction as Fr
from itertools import product, combinations

BASIS = ("1", "r3", "r11", "r33")
MUL = {("1", x): (1, x) for x in BASIS}
MUL.update({(x, "1"): (1, x) for x in BASIS})
MUL.update({("r3", "r3"): (3, "1"), ("r11", "r11"): (11, "1"), ("r33", "r33"): (33, "1"),
            ("r3", "r11"): (1, "r33"), ("r11", "r3"): (1, "r33"),
            ("r3", "r33"): (3, "r11"), ("r33", "r3"): (3, "r11"),
            ("r11", "r33"): (11, "r3"), ("r33", "r11"): (11, "r3")})

def el(**kw):
    return {k: Fr(v) for k, v in kw.items() if v != 0}

def add(a, b, s=1):
    r = dict(a)
    for k, v in b.items():
        r[k] = r.get(k, 0) + s * v
        if r[k] == 0:
            del r[k]
    return r

def mul(a, b):
    r = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            c, k = MUL[(k1, k2)]
            r[k] = r.get(k, 0) + c * v1 * v2
    return {k: v for k, v in r.items() if v != 0}

def scal(c, a):
    return {k: Fr(c) * v for k, v in a.items() if c != 0}

def vadd(P, Q, s=1):
    return (add(P[0], Q[0], s), add(P[1], Q[1], s))

def vscal(c, P):
    return (scal(c, P[0]), scal(c, P[1]))

def d2(P, Q):
    dx, dy = add(P[0], Q[0], -1), add(P[1], Q[1], -1)
    return add(mul(dx, dx), mul(dy, dy))

ZERO = {}
A = (ZERO, ZERO)
D = (el(r11=Fr(1, 2)), el(**{"1": Fr(1, 2)}))
G = (el(r11=Fr(1, 2)), el(**{"1": Fr(-1, 2)}))
# unit normal of D: D/|D| = (sqrt11/(2 sqrt3), 1/(2 sqrt3)); normal (-1/(2 sqrt3), sqrt11/(2 sqrt3)) = (-r3/6, r33/6)
nD = (el(r3=Fr(-1, 6)), el(r33=Fr(1, 6)))
nG = (el(r3=Fr(1, 6)), el(r33=Fr(1, 6)))     # normal of G = (sqrt11/2, -1/2): (1/(2 sqrt3), sqrt11/(2 sqrt3))
B = vadd(vscal(Fr(1, 2), D), vscal(Fr(1, 2), nD))
C = vadd(vscal(Fr(1, 2), D), vscal(Fr(1, 2), nD), -1)
E = vadd(vscal(Fr(1, 2), G), vscal(Fr(1, 2), nG))
F = vadd(vscal(Fr(1, 2), G), vscal(Fr(1, 2), nG), -1)
NAMES = "ABCDEFG"
PTS = [A, B, C, D, E, F, G]
ONE = el(**{"1": 1})
edges = []
for i, j in combinations(range(7), 2):
    dd = d2(PTS[i], PTS[j])
    if dd == ONE:
        edges.append((i, j))
expected = {("A", "B"), ("A", "C"), ("B", "C"), ("B", "D"), ("C", "D"),
            ("A", "E"), ("A", "F"), ("E", "F"), ("E", "G"), ("F", "G"), ("D", "G")}
got = {(NAMES[i], NAMES[j]) for i, j in edges}
assert got == expected, got
assert d2(nD, (ZERO, ZERO)) == ONE and d2(nG, (ZERO, ZERO)) == ONE
print("(1) unit-distance pairs among the 7 points (exact, Q(sqrt3,sqrt11)):", sorted(got), "-> exactly the 11 spindle edges")
for P, nm in zip(PTS, NAMES):
    print("   ", nm, "=", P)

n = 7
adj = [[False] * n for _ in range(n)]
for i, j in edges:
    adj[i][j] = adj[j][i] = True

def colourings(p, q, fix0=True):
    """all maps V -> Z/p with q <= (c(y)-c(x) mod p) <= p-q on edges (backtracking); c(A)=0 if fix0"""
    res = []
    c = [None] * n
    def ok(v, col):
        for u in range(v):
            if adj[u][v]:
                d = (col - c[u]) % p
                if not (q <= d <= p - q):
                    return False
        return True
    def rec(v):
        if v == n:
            res.append(tuple(c))
            return
        for col in ([0] if (fix0 and v == 0) else range(p)):
            if ok(v, col):
                c[v] = col
                rec(v + 1)
        c[v] = None
    rec(0)
    return res

alpha = max(len(S) for r in range(1, n + 1) for S in combinations(range(n), r)
            if all(not adj[i][j] for i, j in combinations(S, 2)))
chi = min(k for k in range(1, 8) if colourings(k, 1))
print(f"(2) alpha = {alpha}, chi = {chi}")
assert alpha == 2 and chi == 4

ok_fracs = []
for p in range(2, 15):
    for q in range(1, p // 2 + 1):
        if 2 * q <= p and Fr(p, q) <= 4:
            if colourings(p, q):
                ok_fracs.append(Fr(p, q))
least = min(ok_fracs)
print(f"(3) least p/q (p <= 14) with a (p,q)-colouring: {least}; fractions < 7/2 with a colouring: "
      f"{sorted(set(f for f in ok_fracs if f < Fr(7, 2)))}")
assert least == Fr(7, 2)

col = (0, 2, 4, 6, 3, 5, 1)
assert all(2 <= (col[j] - col[i]) % 7 <= 5 for i, j in edges)
print("(4) the draft's colouring A..G = 0,2,4,6,3,5,1 is a (7,2)-colouring: OK")

def has_directed_cycle(arcs):
    out = {v: [] for v in range(n)}
    for x, y in arcs:
        out[x].append(y)
    colour = [0] * n
    def dfs(v):
        colour[v] = 1
        for w in out[v]:
            if colour[w] == 1 or (colour[w] == 0 and dfs(w)):
                return True
        colour[v] = 2
        return False
    return any(colour[v] == 0 and dfs(v) for v in range(n))

allc = colourings(7, 2, fix0=False)
bad = 0
for c in allc:
    arcs = []
    for i, j in edges:
        if (c[j] - c[i]) % 7 == 2:
            arcs.append((i, j))
        if (c[i] - c[j]) % 7 == 2:
            arcs.append((j, i))
    if not has_directed_cycle(arcs):
        bad += 1
print(f"(5) number of (7,2)-colourings: {len(allc)}; without a tight cycle: {bad}")
assert bad == 0
print("MOSER: chi_c = 7/2 confirmed (and every (7,2)-colouring has a tight cycle)")
