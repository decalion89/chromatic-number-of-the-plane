"""The local relaxation at any ratio, and why 22/5 differs from 9/2.

Minty at r = p/q asks every cycle for at least |C|q/p edges on its minority
side.  Comparing the two ratios in question:

    length   3  4  5  6  7  8  9  10
    at 9/2   1  1  2  2  2  2  2   3
    at 22/5  1  1  2  2  2  2  3   3

The ONLY difference is at length 9 -- 22/5 demands three where 9/2 demands
two -- and the cycles that hold G at 9/2 are exactly 9-cycles with two.  So
this is the instrument the pending question wants, and it is far cheaper
than encoding 22 positions per vertex.

A relaxation over a SUBSET of cycles is still a relaxation, so it does not
need all the 9-cycles: any collection of them that comes out UNSAT is a
proof.  Enumeration is therefore capped, and the cap costs only strength,
never soundness.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, itertools, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1]
R = Fr(*map(int, sys.argv[2].split("/")))
MAXLEN = int(sys.argv[3]) if len(sys.argv) > 3 else 9
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 400000
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
p, q = R.numerator, R.denominator
need = {L: -(-L * q // p) for L in range(3, MAXLEN + 1)}
print(f"{CAR}: {n} pts, {len(E)} edges; ratio {R} = {float(R):.4f}", flush=True)
print(f"minority side needed by length: {need}  [{time.time()-t0:.0f}s]",
      flush=True)


def lits(cyc):
    out = []
    for x, y in zip(cyc, cyc[1:] + [cyc[0]]):
        v = eid[(min(x, y), max(x, y))]
        out.append(v if x < y else -v)
    return out


def atleast(ls, k):
    return [list(c) for c in itertools.combinations(ls, len(ls) - k + 1)]


s = Solver(name="cd15")
counts, ncl = defaultdict(int), 0
seen = set()
random.seed(2)
PERV = max(1, CAP // n)


def emit(path):
    global ncl
    L = len(path)
    k = need[L]
    if k < 1:
        return
    ls = lits(list(path))
    for c0 in atleast(ls, k) + atleast([-x for x in ls], k):
        s.add_clause(c0)
        ncl += 1
    counts[L] += 1


# Two earlier versions wasted the budget.  The first spent it all on three
# start vertices; the second spread it per vertex but let depth-first search
# eat it with 9-cycles, admitting 241 triangles of 2840 and 666 five-cycles
# of 63604.  Enumerate by LENGTH instead: the short cycles are few and are
# taken whole, and only the long ones are capped.
for L in range(3, MAXLEN + 1):
    if need[L] < 1:
        continue
    cap = CAP if L <= 5 else CAP // 2
    taken = 0
    for start in range(n):
        if taken >= cap:
            break
        stack = [(start, [start])]
        while stack:
            v, path = stack.pop()
            if len(path) == L:
                if start in adj[v]:
                    key = frozenset(path)
                    if key not in seen:
                        seen.add(key)
                        emit(list(path))
                        taken += 1
                        if taken >= cap:
                            break
                continue
            for u in adj[v]:
                if u not in path and u > start:
                    stack.append((u, path + [u]))
    print(f"   length {L}: {counts[L]} cycles, need {need[L]}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"cycles used {dict(sorted(counts.items()))}, {ncl} clauses on "
      f"{len(E)} variables  [{time.time()-t0:.0f}s]", flush=True)
ok = s.solve()
print(("SATISFIABLE -- no conclusion" if ok else
       f"UNSATISFIABLE -- chi_c({CAR}) > {R} = {float(R):.4f} PROVED"),
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
