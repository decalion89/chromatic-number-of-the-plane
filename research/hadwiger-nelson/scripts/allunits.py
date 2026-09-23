"""Every unit vector of an edge module, not just the ones the graph uses.

A unit vector m of the module M has |sigma(m)|^2 = 1 under all 2^r real
embeddings sigma of the coordinate field, so its trace norm
T(m) = sum_sigma |sigma(m)|^2 equals the field degree.  T is a positive
definite form on M: LLL-reduce it, enumerate every lattice point with
T <= degree (Fincke-Pohst), and keep those with |m|^2 = 1 exactly.
"""
import sys, json, math, itertools, time
from fractions import Fraction as Fr
from math import gcd
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
NAME = sys.argv[1]; t0 = time.time()
d = json.load(open(NAME if NAME.startswith("/") else f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
raw = {}
for a, b in g.edges():
    v = tuple((P[b].x - P[a].x).c) + tuple((P[b].y - P[a].y).c); raw[v] = True
den = 1
for v in raw:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(q) * den) for q in v) for v in raw)})
B = echelon(E); r = len(B); D = F.dim
print(f"{NAME}: {len(E)} directions, module rank {r}, field degree {D}", flush=True)
# trace form, exactly: the rational part of x^2 + y^2 is sum_m p_m (x_m^2 + y_m^2)
wt = list(F._prod) * 2
Bint = [list(b) for b in B]
G = [[Fr(sum(wt[k] * Bint[i][k] * Bint[j][k] for k in range(2 * D)), den * den) for j in range(r)] for i in range(r)]
def lll_exact(G, delta=Fr(99, 100)):
    n = len(G); U = [[int(i == j) for j in range(n)] for i in range(n)]
    G = [row[:] for row in G]
    def gso():
        mu = [[Fr(0)] * n for _ in range(n)]; Bn = [Fr(0)] * n
        for i in range(n):
            for j in range(i):
                mu[i][j] = (G[i][j] - sum(mu[j][k] * mu[i][k] * Bn[k] for k in range(j))) / Bn[j]
            Bn[i] = G[i][i] - sum(mu[i][k] ** 2 * Bn[k] for k in range(i))
        return mu, Bn
    def addrow(k, j, q):          # b_k <- b_k - q b_j
        for t in range(n): U[k][t] -= q * U[j][t]
        for t in range(n): G[k][t] -= q * G[j][t]
        for t in range(n): G[t][k] -= q * G[t][j]
    def swap(k):
        U[k], U[k - 1] = U[k - 1], U[k]
        G[k], G[k - 1] = G[k - 1], G[k]
        for row in G: row[k], row[k - 1] = row[k - 1], row[k]
    k = 1
    while k < n:
        mu, Bn = gso()
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                addrow(k, j, q); mu, Bn = gso()
        if Bn[k] >= (delta - mu[k][k - 1] ** 2) * Bn[k - 1]: k += 1
        else: swap(k); k = max(k - 1, 1)
    return G, U
Gr, Ur = lll_exact(G)
Gr = np.array([[float(x) for x in row] for row in Gr])
print(f"  exact LLL done; reduced diagonal {np.round(np.diag(Gr), 4)}   [{time.time()-t0:.0f}s]", flush=True)
# Fincke-Pohst on Gr
n = r; bound = 1 + 1e-9
L = np.linalg.cholesky(Gr)                     # Gr = L L^T
R = L.T                                        # upper triangular: q(z) = |R z|^2
sols = []
z = [0] * n
def rec(i, partial):
    # choose z[i] given z[i+1..]
    s = sum(R[i, j] * z[j] for j in range(i + 1, n))
    rem = bound - partial
    if rem < -1e-9: return
    rad = math.sqrt(max(rem, 0)) / abs(R[i, i])
    c = -s / R[i, i]
    for zi in range(math.ceil(c - rad - 1e-9), math.floor(c + rad + 1e-9) + 1):
        z[i] = zi
        val = (R[i, i] * zi + s) ** 2
        if i == 0:
            if partial + val <= bound and any(z): sols.append(tuple(z))
        else:
            rec(i - 1, partial + val)
    z[i] = 0
rec(n - 1, 0.0)
print(f"  {len(sols)} lattice points with trace norm <= 1 (units have exactly 1)   [{time.time()-t0:.0f}s]", flush=True)
# back to the echelon basis, then ambient vectors, exact unit test
one = F.rational(1)
units = set()
for zz in sols:
    coef = [sum(int(Ur[j][i]) * zz[j] for j in range(n)) for i in range(n)]   # z_reduced -> original basis coefficients
    vec = [sum(coef[i] * Bint[i][k] for i in range(n)) for k in range(2 * D)]
    x = F.element([Fr(vec[m], den) for m in range(D)]); y = F.element([Fr(vec[D + m], den) for m in range(D)])
    if x * x + y * y == one:
        units.add(max(tuple(vec), tuple(-t for t in vec)))
print(f"  unit vectors of the module: {2 * len(units)} ({len(units)} directions) against the graph's {len(E)}   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"source": NAME, "den": den, "directions": [list(u) for u in sorted(units)]},
          open("/home/user/darwin-50/research/hadwiger-nelson/scripts/allunits_" + NAME.split("/")[-1], "w"))
