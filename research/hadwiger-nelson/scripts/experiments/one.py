"""One circular clique, one answer: does this graph map to K(p/q)?

chi_c of a subgraph is a lower bound for chi_c of the whole plane, so a
single UNSAT here is a quantitative statement about R^2 -- sharper than the
"> 4" that chi >= 5 already gives.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1]
R = Fr(*map(int, sys.argv[2].split("/")))
t0 = time.time()
P = (build_G(K, as_graph=False) if CAR == "G"
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
assert b.overflow_headroom(rows) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for v in range(n):
    for j in range(p):
        for j2 in range(j + 1, p):
            cls.append([-(1 + v * p + j), -(1 + v * p + j2)])
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
print(f"{CAR}: {n} pts, {len(E)} edges -> K({p}/{q}) = {float(R):.4f}; "
      f"{n*p} vars, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15", bootstrap_with=cls)
ok = s.solve()
print(f"{'MAPS  -> chi_c <= ' + str(R) if ok else 'DOES NOT MAP  -> chi_c > ' + str(R) + f'  ({float(R):.4f})'}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
