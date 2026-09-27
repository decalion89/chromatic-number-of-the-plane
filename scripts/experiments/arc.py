"""The arc bound: de Grey's palette cap, but measuring spread, not count.

His lemma says a centre and its ring take at most two of four colours -- a
bound on how MANY colours a set can display.  On the circle there is a finer
quantity: not how many positions the set occupies but how wide an arc they
span.  Define A(S) as the largest, over all homomorphisms to K(p/q), of the
shortest arc containing the image of S.  It is a cap in exactly de Grey's
sense -- monotone downward in the carrier, so it travels to supergraphs --
and it says strictly more than a colour count.

The reason to want it: if A(S) < 2q then every two points of S land within
2q of each other on the circle, which is closer than adjacency allows.  So S
is forced INDEPENDENT in every homomorphism.  A forced independent set that
contains two points at distance one is an outright contradiction, with no
rotation argument needed at all.  The colour-counting cap never yields that
directly; it only ever says "few colours", and few is not one.

Computing it: the image of S fits in an arc of length L exactly when the
circle has a run of p-L consecutive positions that S misses.  So ask the
solver for a homomorphism with NO such run -- one clause per starting
position, over auxiliary variables marking the positions S occupies.  UNSAT
at L means every homomorphism confines S to some arc of length L, and the
least such L is A(S).
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
CARRIER = sys.argv[1]
RATIO = Fr(*map(int, sys.argv[2].split("/")))
MODE = sys.argv[3] if len(sys.argv) > 3 else "nbhd"
TOP = int(sys.argv[4]) if len(sys.argv) > 4 else 8
t0 = time.time()
P = (build_G(K, as_graph=False) if CARRIER == "G"
     else pickle.load(open(SC + CARRIER, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
assert b.overflow_headroom(rows) < 1.0
adj = defaultdict(set)
for a, c in fast_edges_complete(b, rows):
    adj[a].add(c)
    adj[c].add(a)
n = len(P)
p, q = RATIO.numerator, RATIO.denominator
print(f"{CARRIER}: {n} points, {sum(len(v) for v in adj.values())//2} edges;"
      f" ratio {RATIO} = {float(RATIO):.4f}   [{time.time()-t0:.0f}s]",
      flush=True)

base = [[1 + v * p + j for j in range(p)] for v in range(n)]
for v in range(n):
    for j in range(p):
        for j2 in range(j + 1, p):
            base.append([-(1 + v * p + j), -(1 + v * p + j2)])
for v in range(n):
    for u in adj[v]:
        if u > v:
            for j in range(p):
                for d in range(-(q - 1), q):
                    base.append([-(1 + v * p + j),
                                 -(1 + u * p + (j + d) % p)])
USED = n * p + 1
print(f"base: {len(base)} clauses   [{time.time()-t0:.0f}s]", flush=True)


def arc_bound(S):
    """Least L such that every homomorphism puts S inside an arc of L."""
    S = list(S)
    anchor = [[-(USED + j)] + [1 + s * p + j for s in S] for j in range(p)]
    for L in range(1, p + 1):
        gap = p - L
        if gap <= 0:
            return p
        cls = base + anchor + [[USED + (i + t) % p for t in range(gap)]
                               for i in range(p)]
        s = Solver(name="cd15", bootstrap_with=cls)
        ok = s.solve()
        s.delete()
        if not ok:
            return L
    return p


deg = sorted(range(n), key=lambda v: -len(adj[v]))
if MODE == "nbhd":
    targets = [("v%d+N(v) deg %d" % (v, len(adj[v])), [v] + list(adj[v]))
               for v in deg[:TOP]]
elif MODE == "ring":
    targets = [("N(v%d) deg %d" % (v, len(adj[v])), list(adj[v]))
               for v in deg[:TOP]]
else:
    targets = [("whole carrier", list(range(n)))]
print(f"{len(targets)} targets; forced independent when A < 2q = {2*q}",
      flush=True)
for name, S in targets:
    A = arc_bound(S)
    pairs = sum(1 for i, a in enumerate(S) for c in S[i + 1:] if c in adj[a])
    verdict = ("FORCED INDEPENDENT" if A < 2 * q else "free")
    print(f"  |S|={len(S):4d}  A(S) = {A:2d} of {p}   {verdict}"
          f"   ({pairs} unit pairs inside S)   {name}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if A < 2 * q and pairs:
        print(f"  *** CONTRADICTION: S is forced independent but holds "
              f"{pairs} unit pairs ***", flush=True)
print("DONE", flush=True)
