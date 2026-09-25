"""The same question on a big carrier, built clause by clause.

chi_c(R^2) is at least chi_c of any unit-distance graph, and a bigger graph
can only give a bigger bound.  The universe here has 39144 points, so the
clause list would be millions long: it is fed to the solver one clause at a
time and never held in Python, which is what makes the size affordable.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, gc
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR, R = sys.argv[1], Fr(*map(int, sys.argv[2].split("/")))
LIM = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
t0 = time.time()
P = pickle.load(open(SC + CAR, "rb"))[:LIM]
b = IntBasis.covering(P)
rows = b.rows(P)
assert b.overflow_headroom(rows) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
del rows
gc.collect()
n, p, q = len(P), R.numerator, R.denominator
print(f"{CAR}: {n} pts, {len(E)} edges -> K({p}/{q}) = {float(R):.4f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
s = Solver(name="cd15")
for v in range(n):
    s.add_clause([1 + v * p + j for j in range(p)])
    for j in range(p):
        for j2 in range(j + 1, p):
            s.add_clause([-(1 + v * p + j), -(1 + v * p + j2)])
print(f"vertex clauses in  [{time.time()-t0:.0f}s]", flush=True)
for k, (a, c) in enumerate(E):
    for j in range(p):
        for d in range(-(q - 1), q):
            s.add_clause([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
    if k % 50000 == 0:
        print(f"   {k}/{len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)
print(f"all clauses in, solving  [{time.time()-t0:.0f}s]", flush=True)
ok = s.solve()
print(("MAPS  -> chi_c <= " + str(R)) if ok else
      ("DOES NOT MAP  -> chi_c(R^2) > " + str(R) + f" = {float(R):.4f}"),
      f" [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
