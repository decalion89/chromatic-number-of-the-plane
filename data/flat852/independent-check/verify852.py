"""Independent check of core852a.json (written apart from the growth agent's code): exact unit distances in
Q(zeta21) via sympy's cyclotomic polynomial, 50-digit numerics, distinct vertices, a real triangle, and my own
4-colouring CNF (variable c*n + v + 1)."""
import json, sys, itertools
import sympy as sp
import mpmath as mp
mp.mp.dps = 50
g = json.load(open(sys.argv[1]))
V = [tuple(v["exact7"]) for v in g["vertices"]]; n = len(V)
E = [(a, b) for a, b, j in g["edges"]]
x = sp.symbols("x"); PHI = sp.Poly(sp.cyclotomic_poly(21, x), x)
assert PHI.degree() == 12
def poly(c): return sp.Poly(sum(int(ci) * x**i for i, ci in enumerate(c)), x)
def conj(p):  # z -> z^20 = z^-1 on Q(zeta21), then reduce
    return sp.Poly(sum(p.coeff_monomial(x**i) * x**((21 - i) % 21) for i in range(p.degree() + 1)), x).rem(PHI)
assert len(set(V)) == n, "repeated vertex"
bad = 0
for a, b in E:
    d = poly([vb - va for va, vb in zip(V[a], V[b])])      # 7 * (v_b - v_a)
    if d.mul(conj(d)).rem(PHI) != sp.Poly(49, x): bad += 1
print("edges", len(E), "exact unit-distance failures", bad)
z = mp.exp(2j * mp.pi / 21)
P = [sum(c * z**i for i, c in enumerate(v)) / 7 for v in V]
worst = max(abs(abs(P[b] - P[a]) - 1) for a, b in E)
print("50-digit max | |d|-1 | over edges:", mp.nstr(worst, 5))
mind = min(abs(P[i] - P[k]) for i in range(n) for k in range(i + 1, n))
print("min distance between vertices:", mp.nstr(mind, 5))
adj = set(map(frozenset, E)); t = g["triangle_fixed"]
assert all(frozenset(p) in adj for p in itertools.combinations(t, 2)), "not a triangle"
col = json.load(open(sys.argv[2]))
if isinstance(col, dict):
    col = col.get("colouring", col.get("colours", col))
if isinstance(col, dict):
    col = [col[str(i)] if str(i) in col else col[i] for i in range(n)]
print("5-colouring proper:", all(col[a] != col[b] for a, b in E), "colours used:", len(set(col)))
var = lambda v, c: c * n + v + 1
cl = [[var(v, c) for c in range(4)] for v in range(n)]
cl += [[-var(a, c), -var(b, c)] for a, b in E for c in range(4)]
cl += [[var(t[i], i)] for i in range(3)]
with open("mine852.cnf", "w") as f:
    f.write(f"c independent 4-colouring CNF of core852a (variable c*n+v+1), triangle {t} fixed to 0,1,2\n")
    f.write(f"p cnf {4 * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
print("wrote mine852.cnf", 4 * n, "vars", len(cl), "clauses; bad edges must be 0 to use it")
