"""Forced-monochromatic independent sets: the strongest lemma shape there is.

De Grey's cap is 2 of 4 on a seven-point set, and densifying Sa eightfold to
3501 points leaves it at exactly 2.  It cannot do better: the set contains
unit edges, so a cap of 1 would put adjacent points in one colour.  His
lemma is already optimal FOR ITS SHAPE, and no carrier can improve it.

The only stronger lemma is a cap of 1 on an INDEPENDENT set -- every point
forced to the same colour in every proper colouring.  That hands over the
forced pair directly, where the bite needs a disjunction and a case analysis
to extract one, and a forced-monochromatic set of three or more would be
something the construction has never had.

Forced pairs are known to exist here; this asks for bigger ones.  The test
is one solver call: assert two colours appear on S and look for UNSAT.
"""
import sys, time, pickle, itertools
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "Sa"
KC = int(sys.argv[2]) if len(sys.argv) > 2 else 4
SRC = int(sys.argv[3]) if len(sys.argv) > 3 else 200
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve(), f"{CAR} is not {KC}-colourable"
print(f"{CAR}: {n} pts, {len(E)} edges, k = {KC}  [{time.time()-t0:.0f}s]",
      flush=True)
# The first version built KC^2 selector variables per pair, which is
# hundreds of thousands of them across a scan.  Colour permutation makes it
# unnecessary: if u and v CAN differ then they can differ with u at colour 0
# and v at colour 1, so the whole question is one solve under two
# assumptions and no added clauses at all.
def forced_same(u, v):
    return not s.solve(assumptions=[1 + u * KC, 1 + v * KC + 1])


# And a forced-monochromatic set is just a clique of forced pairs, since
# "all of S one colour" is "every pair of S agrees".  So triples and larger
# come out of the pair relation for free.
pairs = 0
found = []
cands = [v for v in range(min(n, SRC))]
print(f"scanning pairs from {len(cands)} source vertices"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for i, u in enumerate(cands):
    for v in range(u + 1, n):
        if v in adj[u]:
            continue
        pairs += 1
        if forced_same(u, v):
            found.append((u, v))
            print(f"   forced pair {u},{v}  [{time.time()-t0:.0f}s]",
                  flush=True)
    if i % 20 == 19:
        print(f"   {i+1} sources, {pairs} pairs, {len(found)} forced"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if len(found) >= 40 or time.time() - t0 > 900:
        break
print(f"\n{len(found)} forced pairs from {pairs} tested"
      f"  [{time.time()-t0:.0f}s]", flush=True)
# Larger forced-monochromatic sets are cliques in the forced-pair relation.
if found:
    nb = defaultdict(set)
    for u, v in found:
        nb[u].add(v)
        nb[v].add(u)
    triples = []
    for u, v in found:
        for w in nb[u] & nb[v]:
            t = tuple(sorted((u, v, w)))
            if t not in triples:
                triples.append(t)
    print(f"forced-monochromatic triples among them: {len(triples)}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for t in triples[:6]:
        print(f"   *** triple {t} ***", flush=True)
    if not triples:
        print("   none -- the forced pairs do not chain, so the strongest "
              "lemma available here is still a pair", flush=True)
print("DONE", flush=True)
