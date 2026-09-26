"""Rank the graphs this project has built by how close they are to six.

Until now every one of them answered "5" and there was no way to tell them
apart.  The circular chromatic number separates them: it is monotone under
adding points, it lives in (4, 5] for all of them, and chi_c > 5 is exactly
the goal.  So the graph with the highest chi_c is the best place to grow
from, and that question has never been askable before.

Full bisection per graph is dear, so this ranks by a single discriminating
question: does the graph map to K(p/q) for one fixed ratio?  A graph that
fails where another succeeds has the strictly higher chi_c.  Sweep a few
ratios and the ranking comes out of the pattern of failures.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, glob, os
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
RATIOS = [Fr(*map(int, r.split("/"))) for r in sys.argv[1].split(",")]
MAXN = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
t0 = time.time()


def edges_of(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        return None
    return sorted(set((min(a, c), max(a, c))
                      for a, c in fast_edges_complete(b, r)))


def maps_to(n, E, p, q):
    cls = [[1 + v * p + i for i in range(p)] for v in range(n)]
    for v in range(n):
        for i in range(p):
            for j in range(i + 1, p):
                cls.append([-(1 + v * p + i), -(1 + v * p + j)])
    for a, c in E:
        for i in range(p):
            for d in range(-(q - 1), q):
                cls.append([-(1 + a * p + i), -(1 + c * p + (i + d) % p)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    return ok


cands = [("G", build_G(K, as_graph=False))]
for f in sorted(glob.glob(SC + "*.pkl")):
    try:
        P = pickle.load(open(f, "rb"))
    except Exception:
        continue
    if isinstance(P, list) and 100 <= len(P) <= MAXN:
        cands.append((os.path.basename(f), P))
print(f"{len(cands)} graphs; ratios "
      f"{', '.join(str(r) for r in RATIOS)}", flush=True)
rows = []
for name, P in cands:
    E = edges_of(P)
    if E is None or not E:
        print(f"  {name:24s} {len(P):6d} pts  -- skipped", flush=True)
        continue
    res, lowest = [], None
    for r in RATIOS:
        p, q = r.numerator, r.denominator
        if p * len(P) > 400000:
            res.append("-")
            continue
        ok = maps_to(len(P), E, p, q)
        res.append("m" if ok else "X")
        if ok and lowest is None:
            lowest = r
    rows.append((lowest if lowest else Fr(99), name, len(P), len(E),
                 "".join(res)))
    print(f"  {name:24s} {len(P):6d} pts {len(E):7d} edges  "
          f"{''.join(res)}   chi_c <= {lowest if lowest else '>all'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("\n=== ranked, hardest first ===", flush=True)
for lo, name, n, m, res in sorted(rows, reverse=True):
    print(f"  chi_c > {RATIOS[0] if res[0]=='X' else '4'}"
          f"   <= {lo if lo != Fr(99) else 'ALL FAIL'}   {name}"
          f"  ({n} pts, {m} edges)", flush=True)
print("DONE", flush=True)
