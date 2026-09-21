"""Does ANY set of m points fall below the generic palette?  Decided, not
searched.

Hill climbing over sets failed its own calibration, and the reason is
structural: a tight set is an isolated minimum, since swapping any one of its
points out restores the full palette, so a local search has no approach path.
The question needs an exact method.

It has one.  A set S fails to be tight as soon as ONE colouring shows all k
colours on it, so a finite sample of colourings gives a necessary condition
that is easy to encode: pick S with a variable per vertex, and for every
sampled colouring require that some colour is missing from S.  Solving that
returns a candidate.  Verify the candidate exactly; if it really is tight, the
search is over, and if it is not, the verification hands back the very
colouring that breaks it, which goes into the sample and rules out that
candidate and many others at once.

The loop is a decision procedure, not a heuristic.  A truly tight set
satisfies the condition for every sample, so when the encoding becomes
unsatisfiable no tight set of that size exists at all -- an absence proof
rather than a failure to find.

Calibrated at four colours, where the answer is known to be inside Sa.
"""
import sys, time, random, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "Sa"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
M = int(sys.argv[3]) if len(sys.argv) > 3 else 7
TARGET = int(sys.argv[4]) if len(sys.argv) > 4 else k - 1
ROUNDS = int(sys.argv[5]) if len(sys.argv) > 5 else 4000
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
if SEED.endswith(".pkl"):
    P = pickle.load(open(SC + SEED, "rb"))
else:
    P = {"Sa": build_Sa, "Y": build_Y,
         "G": lambda f: build_G(f, as_graph=False)}[SEED](K)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{SEED}: {n} points, {len(E)} edges, k={k}; looking for {M} points "
      f"with palette <= {TARGET}  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve(), f"{SEED} is not {k}-colourable"
rng = random.Random(2654435761)


def a_colouring():
    col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                           for w in range(n * k)])
    col_solver.solve()
    m = col_solver.get_model()
    return [next(c for c in range(k) if m[v * k + c] > 0) for v in range(n)]


def witness(S):
    """A colouring showing more than TARGET colours on S, or None."""
    pool = IDPool(start_from=n * k + 1)
    ind = [pool.id(("u", c)) for c in range(k)]
    ext = list(cls)
    for c in range(k):
        ext.append([-ind[c]] + [1 + v * k + c for v in S])
    ext += list(CardEnc.atleast(lits=ind, bound=TARGET + 1, vpool=pool,
                                encoding=EncType.seqcounter))
    s = Solver(name="cd15", bootstrap_with=ext)
    ok = s.solve()
    out = None
    if ok:
        m = s.get_model()
        out = [next(c for c in range(k) if m[v * k + c] > 0)
               for v in range(n)]
    s.delete()
    return out


samples = [a_colouring() for _ in range(40)]
print(f"seeded with {len(samples)} colourings  [{time.time()-t0:.0f}s]",
      flush=True)
for rnd in range(ROUNDS):
    pool = IDPool(start_from=n + 1)
    y = list(range(1, n + 1))                 # y[v-1] = vertex v-1 is in S
    f = list(CardEnc.equals(lits=y, bound=M, vpool=pool,
                            encoding=EncType.seqcounter))
    for si, c in enumerate(samples):
        byc = defaultdict(list)
        for v in range(n):
            byc[c[v]].append(v)
        miss = [pool.id(("z", si, q)) for q in range(k)]
        # at most TARGET colours present means at least k-TARGET are missing
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
        print(f"\nUNSATISFIABLE after {len(samples)} colourings: NO set of "
              f"{M} points in {SEED} has palette <= {TARGET} at k={k}.  "
              f"This is an absence proof, not a failed search."
              f"  [{time.time()-t0:.0f}s]", flush=True)
        print("DONE", flush=True)
        raise SystemExit(0)
    S = [v for v in range(n) if mod[v] > 0]
    w = witness(S)
    if w is None:
        p = ring_palette_bound(cls, n * k, S, k)
        print(f"\n*** TIGHT SET FOUND: {M} points, palette {p} <= {TARGET} "
              f"at k={k} after {len(samples)} colourings: {S}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        with open(SC + f"tightset_{SEED}_{k}_{M}.pkl", "wb") as fh:
            pickle.dump([P[v] for v in S], fh)
        print("DONE", flush=True)
        raise SystemExit(0)
    samples.append(w)
    if rnd % 50 == 0:
        print(f"   round {rnd}: {len(samples)} colourings"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nran out of rounds with {len(samples)} colourings -- undecided"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
