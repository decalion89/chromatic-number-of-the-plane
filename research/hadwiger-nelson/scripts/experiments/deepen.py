"""Look for palette 3 inside the palette-4 set, where the search is small.

The capped set S has palette 4 of 5, and everything now turns on whether a
palette of 3 is reachable.  With 4 the pigeonhole forces a class of |T|/4,
and the largest set avoiding distances 1 and sqrt3 in a lattice-like carrier
is the spacing-2 sublattice at density exactly 1/4 -- level with the bar, so
one distance can never block and the disjunction needs two.  With 3 the bar
rises to |T|/3 while that sublattice stays at 1/4, and sqrt3 alone blocks.
The whole collapse hangs on one step of depth.

Looking for it in the whole graph is expensive.  Looking for it inside S is
not: any subset of S already has palette at most 4, so the question is only
whether some subset drops to 3, and that is a search over sixty points rather
than fifteen hundred.  The same decision procedure applies, and
unsatisfiability here is a proof that no subset of S reaches depth 3.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
M = int(sys.argv[2]) if len(sys.argv) > 2 else 6
TARGET = int(sys.argv[3]) if len(sys.argv) > 3 else 3
SEEDN = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
SCOPE = sys.argv[5] if len(sys.argv) > 5 else "S"
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
S, _ = pickle.load(open(SC + f"grown_{k}_4.pkl", "rb"))
S = sorted(S)
if SCOPE == "S":
    scope = S
else:                                   # S plus everything one hop out
    ext = set(S)
    for v in S:
        ext |= adj[v]
    scope = sorted(ext)
print(f"scope: {len(scope)} points (S has {len(S)}), looking for {M} of them "
      f"at palette <= {TARGET} of {k}  [{time.time()-t0:.0f}s]", flush=True)
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve()
rng = random.Random(6700417)


def a_colouring():
    col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                           for w in range(n * k)])
    col_solver.solve()
    mo = col_solver.get_model()
    return [next(c for c in range(k) if mo[v * k + c] > 0) for v in range(n)]


def witness(T):
    pool = IDPool(start_from=n * k + 1)
    ind = [pool.id(("u", c)) for c in range(k)]
    ext = list(cls)
    for c in range(k):
        ext.append([-ind[c]] + [1 + v * k + c for v in T])
    ext += list(CardEnc.atleast(lits=ind, bound=TARGET + 1, vpool=pool,
                                encoding=EncType.seqcounter))
    s = Solver(name="cd15", bootstrap_with=ext)
    ok = s.solve()
    out = None
    if ok:
        mo = s.get_model()
        out = [next(c for c in range(k) if mo[v * k + c] > 0)
               for v in range(n)]
    s.delete()
    return out


samples = [a_colouring() for _ in range(SEEDN)]
print(f"seeded with {len(samples)} colourings  [{time.time()-t0:.0f}s]",
      flush=True)
for rnd in range(30000):
    pool = IDPool(start_from=n + 1)
    y = [v + 1 for v in scope]
    f = list(CardEnc.equals(lits=y, bound=M, vpool=pool,
                            encoding=EncType.seqcounter))
    for si, c in enumerate(samples):
        byc = defaultdict(list)
        for v in scope:
            byc[c[v]].append(v)
        miss = [pool.id(("z", si, q)) for q in range(k)]
        f += list(CardEnc.atleast(lits=miss, bound=k - TARGET, vpool=pool,
                                  encoding=EncType.seqcounter))
        for q in range(k):
            for v in byc[q]:
                f.append([-miss[q], -(v + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    mod = s.get_model() if ok else None
    s.delete()
    if not ok:
        print(f"\nUNSATISFIABLE after {len(samples)} colourings: no {M} "
              f"points of this scope reach palette <= {TARGET}.  An absence "
              f"proof.  [{time.time()-t0:.0f}s]", flush=True)
        break
    T = [v for v in scope if mod[v] > 0]
    w = witness(T)
    if w is None:
        p = ring_palette_bound(cls, n * k, T, k)
        print(f"\n*** DEPTH {p} FOUND on {M} points: {T}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        with open(SC + f"deep_{k}_{TARGET}_{M}.pkl", "wb") as fh:
            pickle.dump((T, [P[v] for v in T]), fh)
        break
    samples.append(w)
    if rnd % 200 == 0 and rnd:
        print(f"   round {rnd}: {len(samples)} colourings"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
