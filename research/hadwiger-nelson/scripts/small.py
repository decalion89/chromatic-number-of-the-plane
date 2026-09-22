"""chi_c of the small unit-distance graphs, and a unit test of the encoder.

Odd cycles have chi_c(C_{2k+1}) = 2 + 1/k exactly, and complete graphs have
chi_c(K_m) = m, so running the same code on those says whether it is
measuring what it claims to measure before any of its answers about the
plane are believed.
"""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.geometry import DEGREY_FIELD as K, Point
from hn.graph import build_graph
from pysat.solvers import Solver

half, r3, r11 = K.rational(Fr(1, 2)), K.sqrt(3), K.sqrt(11)
A = Point(K.zero(), K.zero())
D = Point(r3, K.zero())
B = Point(r3 * half, half)
C = Point(r3 * half, -half)
c, s = K.rational(Fr(5, 6)), r11 * K.rational(Fr(1, 6))
sp = [A, B, C, D] + [Point(p.x * c - p.y * s, p.x * s + p.y * c)
                     for p in (B, C, D)]


def maps(n, E, p, q):
    cl = [[1 + v * p + j for j in range(p)] for v in range(n)]
    for a, b in E:
        for j in range(p):
            for d in range(-(q - 1), q):
                cl.append([-(1 + a * p + j), -(1 + b * p + (j + d) % p)])
    cl += [[-(1 + j)] for j in range(1, p)] + [[1]]
    sv = Solver(name="cd15", bootstrap_with=cl)
    ok = sv.solve()
    sv.delete()
    return ok


def chi_c(name, n, E, hi=8):
    cands = sorted({Fr(p, q) for q in range(1, 13) for p in range(2, hi * q + 1)
                    if Fr(p, q).denominator == q and 2 <= Fr(p, q) <= hi})
    for r in cands:
        if maps(n, E, r.numerator, r.denominator):
            print(f"  {name:22s} n={n:3d} m={len(E):3d}   chi_c = {r}"
                  f" = {float(r):.4f}")
            return r
    print(f"  {name:22s} maps to nothing up to {hi}")
    return None


print("unit tests (known values):")
for k in (1, 2, 3, 4):
    m = 2 * k + 1
    chi_c(f"C_{m}  (want 2+1/{k})", m, [(i, (i + 1) % m) for i in range(m)])
for m in (3, 4, 5):
    chi_c(f"K_{m}  (want {m})", m,
          [(i, j) for i in range(m) for j in range(i + 1, m)])
print("\nunit-distance graphs:")
g = build_graph(sp)
E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
chi_c("Moser spindle", g.n, E)
rh = build_graph([A, B, C, D])
chi_c("unit rhombus", rh.n,
      sorted(set((min(a, b), max(a, b)) for a, b in rh.edges())))
