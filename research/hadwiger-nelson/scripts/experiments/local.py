"""A local relaxation of Minty: short cycles only, and a valid lower bound.

chi_c(H) <= 9/2 requires an orientation in which EVERY cycle has at least
2|C|/9 edges on its minority side.  Imposing that over all cycles is
hopeless -- there are exponentially many -- but imposing it over the short
ones is a relaxation, and a relaxation that is UNSAT proves the original is
too.  So:

    triangles and 4-cycles   at least 1 on the minority side (acyclicity)
    5- to 9-cycles           at least 2

UNSAT  =>  chi_c(H) > 9/2, proved.
SAT    =>  nothing; the witness may fail on a longer cycle.

The encoding is tiny.  One boolean per EDGE, not one per vertex per
position, and "at least 2 of 5 true" is five clauses of width 4 -- no
cardinality machinery.  For G that is 7877 variables where the circular
clique encoding wanted 14229, and about 92000 clauses where it wanted
636904.  Since a graph that maps to K(9/2) hands over a satisfying
orientation directly, this can only ever be SAT for such a graph; its use is
on candidates, where an UNSAT is a proof and arrives cheaply.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, itertools
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1]
MAXLEN = int(sys.argv[2]) if len(sys.argv) > 2 else 5
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
eid = {e: i + 1 for i, e in enumerate(E)}
adj = [set() for _ in range(n)]
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
print(f"{CAR}: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]",
      flush=True)


def lits(cyc):
    """Literals true when the edge runs forward along the traversal."""
    out = []
    for x, y in zip(cyc, cyc[1:] + [cyc[0]]):
        v = eid[(min(x, y), max(x, y))]
        out.append(v if x < y else -v)
    return out


def atleast(ls, k):
    """At least k of ls true: no (len-k+1)-subset is entirely false."""
    return [list(c) for c in itertools.combinations(ls, len(ls) - k + 1)]


cls, counts = [], defaultdict(int)
seen = set()
for start in range(n):
    stack = [(start, [start])]
    while stack:
        v, path = stack.pop()
        if len(path) > MAXLEN:
            continue
        for u in adj[v]:
            if u == start and len(path) >= 3:
                key = frozenset(path)
                if len(key) == len(path) and key not in seen:
                    seen.add(key)
                    L = len(path)
                    need = 1 if L <= 4 else 2
                    ls = lits(list(path))
                    cls += atleast(ls, need)
                    cls += atleast([-x for x in ls], need)
                    counts[L] += 1
            elif u not in path and u > start and len(path) < MAXLEN:
                stack.append((u, path + [u]))
print(f"cycles by length {dict(sorted(counts.items()))}; "
      f"{len(cls)} clauses on {len(E)} variables  [{time.time()-t0:.0f}s]",
      flush=True)
s = Solver(name="cd15", bootstrap_with=cls)
ok = s.solve()
print(("SATISFIABLE -- no conclusion (a longer cycle may still fail)" if ok
       else "UNSATISFIABLE -- chi_c > 9/2 PROVED by short cycles alone"),
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
