"""Own check of the remark 'the 76-vertex graph of chi(Q(sqrt11)^2) = 4 has chi_c = 16/5'.
(1) exact check that the stored edges are exactly the unit-distance pairs (points ((a+b r)/D, (c+e r)/D), r^2 = 11);
(2) chi = 4 (SAT for 4, UNSAT for 3);  (3) a (16,5)-colouring exists (verified);  (4) no (P,Q)-colouring for the
largest P/Q < 16/5 with P <= 76 (then chi_c = 16/5 by the theorem that chi_c(H) = p/q with p <= |V(H)|).
SAT: pysat / CaDiCaL, own direct encoding; UNSAT answers are not certified (no proof checked)."""
import os
import json, sys, time
from fractions import Fraction as Fr
from itertools import combinations
from pysat.solvers import Solver
g = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "quadratic_planes", "q11.json")))
D, pts, E = g["D"], g["points"], [tuple(e) for e in g["edges"]]
n = len(pts)
def sqd(P, Q):
    # coordinates (a + b r)/D: difference squared, as (rational part, r part)
    dx = (P[0] - Q[0], P[1] - Q[1]); dy = (P[2] - Q[2], P[3] - Q[3])
    rat = dx[0] ** 2 + 11 * dx[1] ** 2 + dy[0] ** 2 + 11 * dy[1] ** 2
    irr = 2 * dx[0] * dx[1] + 2 * dy[0] * dy[1]
    return (Fr(rat, D * D), Fr(irr, D * D))
unit = sorted(tuple(sorted(e)) for e in combinations(range(n), 2) if sqd(pts[e[0]], pts[e[1]]) == (1, 0))
assert unit == sorted(tuple(sorted(e)) for e in E), "edge list differs from the unit-distance pairs"
print(f"(1) {n} points, {len(E)} edges = exactly the unit-distance pairs (exact arithmetic in Q(sqrt11))")
def hom(p, q, tlim=None):
    var = lambda v, a: v * p + a + 1
    s = Solver(name="cadical153")
    for v in range(n):
        s.add_clause([var(v, a) for a in range(p)])
    bad = [d for d in range(p) if not (q <= d <= p - q)]
    for (u, v) in E:
        for a in range(p):
            for d in bad:
                s.add_clause([-var(u, a), -var(v, (a + d) % p)])
    s.add_clause([var(E[0][0], 0)])
    t = time.time(); ok = s.solve(); dt = time.time() - t
    col = None
    if ok:
        m = set(l for l in s.get_model() if l > 0)
        col = [next(a for a in range(p) if var(v, a) in m) for v in range(n)]
        assert all(q <= (col[v] - col[u]) % p <= p - q for u, v in E)
    s.delete()
    return ok, dt
for (p, q) in [(3, 1), (4, 1), (16, 5)]:
    ok, dt = hom(p, q)
    print(f"K_{p}/{q}: {'SAT (colouring verified)' if ok else 'UNSAT'} [{dt:.1f}s]", flush=True)
cands = sorted({Fr(P, Q) for P in range(2, n + 1) for Q in range(1, P // 2 + 1) if Fr(P, Q) < Fr(16, 5)})
pred = cands[-1]
print("largest P/Q < 16/5 with P <=", n, ":", pred, flush=True)
ok, dt = hom(pred.numerator, pred.denominator)
print(f"K_{pred.numerator}/{pred.denominator}: {'SAT' if ok else 'UNSAT'} [{dt:.1f}s]")
print("=> chi_c = 16/5 (modulo the uncertified UNSAT)" if not ok else "=> chi_c < 16/5 !")
