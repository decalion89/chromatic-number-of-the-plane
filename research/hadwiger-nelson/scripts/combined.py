"""Ask for both halves at once: capped, and 5-chromatic for two distances.

The reduction is settled and both halves exist separately.  A set capped at
palette 4 splits into four classes each independent for distance 1; if its
{1, d} graph needed five colours those classes could not avoid d either, and
some class would hold a monochromatic pair at distance d exactly.  Three
hundred cores have chi 5 and palette 5; the capped sets have palette 4 and
chi 3.  Nothing has had both.

Both conditions are of the same shape -- "for every X, the chosen set does
something" -- so both are CEGAR-able against the same variables.

The cap needs sampled 5-colourings of the carrier: the set must miss a colour
in each.  The chromatic condition needs sampled 4-colourings, and here is the
part that makes it work: the sampled colouring does NOT have to be proper.  If
the chosen set's two-distance graph is not 4-colourable, then no assignment of
four colours to it is proper, so EVERY four-colouring of the carrier restricted
to it leaves a monochromatic edge.  Arbitrary assignments are therefore valid
constraints, and there are as many of them as one cares to draw.

Unsatisfiable is then a proof that no set of that size in that scope has both
properties.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pairsat import pairs_at
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
M = int(sys.argv[2]) if len(sys.argv) > 2 else 13
DSQ = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(3)
CARRIER = sys.argv[4] if len(sys.argv) > 4 else "G"
NCOL = int(sys.argv[5]) if len(sys.argv) > 5 else 400
NFOUR = int(sys.argv[6]) if len(sys.argv) > 6 else 400
HOPS = int(sys.argv[7]) if len(sys.argv) > 7 else 0
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = (build_G(K, as_graph=False) if CARRIER == "G"
     else pickle.load(open(SC + CARRIER, "rb")))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E1 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
D2 = b.D * b.D
s2 = pairs_at(b, r, int(DSQ * D2))
two = sorted(set(E1) | set(s2))
print(f"{CARRIER}: {n} points, {len(E1)} unit + {len(s2)} at distance^2 "
      f"{DSQ}; asking for {M} points capped at {k-1} whose two-distance "
      f"graph needs {k} colours  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E1:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve()
rng = random.Random(3010349)


def a_colouring():
    col_solver.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                           for w in range(n * k)])
    col_solver.solve()
    mo = col_solver.get_model()
    return [next(c for c in range(k) if mo[v * k + c] > 0) for v in range(n)]


def witness_cap(T):
    pool = IDPool(start_from=n * k + 1)
    ind = [pool.id(("u", c)) for c in range(k)]
    ext = list(cls)
    for c in range(k):
        ext.append([-ind[c]] + [1 + v * k + c for v in T])
    ext += list(CardEnc.atleast(lits=ind, bound=k, vpool=pool,
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


def four_colour(T):
    """A proper 4-colouring of T's two-distance graph, or None."""
    loc = {v: i for i, v in enumerate(T)}
    ts = set(T)
    c2 = [[1 + i * 4 + c for c in range(4)] for i in range(len(T))]
    for a, c in two:
        if a in ts and c in ts:
            for col in range(4):
                c2.append([-(1 + loc[a] * 4 + col), -(1 + loc[c] * 4 + col)])
    s = Solver(name="cd15", bootstrap_with=c2)
    ok = s.solve()
    out = None
    if ok:
        mo = s.get_model()
        out = {v: next(c for c in range(4) if mo[loc[v] * 4 + c] > 0)
               for v in T}
    s.delete()
    return out


# A local scope is what makes the encoding affordable: each sampled
# four-assignment needs an indicator per monochromatic pair inside the chosen
# region, and over the whole graph that is a million variables before the
# solver has seen a single clause.
if HOPS:
    adjs = defaultdict(set)
    for a, c in E1:
        adjs[a].add(c)
        adjs[c].add(a)
    seed = pickle.load(open(SC + "core13.pkl", "rb"))[0]
    cur = set(seed)
    for _ in range(HOPS):
        nxt = set(cur)
        for v in cur:
            nxt |= adjs[v]
        cur = nxt
    scope = sorted(cur)
else:
    scope = list(range(n))
ins = set(scope)
two = [(a, c) for a, c in two if a in ins and c in ins]
print(f"scope: {len(scope)} points, {len(two)} two-distance pairs inside",
      flush=True)
fives = [a_colouring() for _ in range(NCOL)]
fours = [[rng.randrange(4) for _ in range(n)] for _ in range(NFOUR)]
print(f"seeded with {len(fives)} proper 5-colourings and {len(fours)} "
      f"arbitrary 4-assignments  [{time.time()-t0:.0f}s]", flush=True)
for rnd in range(30000):
    if rnd % 20 == 0:
        print(f"   round {rnd}: {len(fives)} colourings, {len(fours)} "
              f"assignments  [{time.time()-t0:.0f}s]", flush=True)
    pool = IDPool(start_from=n + 1)
    f = list(CardEnc.equals(lits=[v + 1 for v in scope], bound=M, vpool=pool,
                            encoding=EncType.seqcounter))
    for si, c in enumerate(fives):
        byc = defaultdict(list)
        for v in scope:
            byc[c[v]].append(v)
        miss = [pool.id(("z", si, q)) for q in range(k)]
        f.append(miss)
        for q in range(k):
            for v in byc[q]:
                f.append([-miss[q], -(v + 1)])
    for fi, phi in enumerate(fours):
        lits = []
        for a, c in two:
            if phi[a] == phi[c]:
                e = pool.id(("m", fi, a, c))
                f.append([-e, a + 1])
                f.append([-e, c + 1])
                lits.append(e)
        f.append(lits if lits else [])
        if not lits:
            print("   an assignment makes the demand impossible -- skipped",
                  flush=True)
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    mod = s.get_model() if ok else None
    s.delete()
    if not ok:
        print(f"\nUNSATISFIABLE after {len(fives)} colourings and "
              f"{len(fours)} assignments: no {M} points of {CARRIER} are both "
              f"capped and {k}-chromatic for {{1, sqrt{DSQ}}}.  An absence "
              f"proof.  [{time.time()-t0:.0f}s]", flush=True)
        break
    T = [v for v in scope if mod[v] > 0]
    w = witness_cap(T)
    if w is not None:
        fives.append(w)
        continue
    phi = four_colour(T)
    if phi is not None:
        full = [rng.randrange(4) for _ in range(n)]
        for v, c in phi.items():
            full[v] = c
        fours.append(full)
        continue
    pal = ring_palette_bound(cls, n * k, T, k)
    print(f"\n*** BOTH: {M} points, palette {pal} of {k}, and its "
          f"{{1, sqrt{DSQ}}} graph is not 4-colourable.  In EVERY "
          f"{k}-colouring of {CARRIER} two of them at distance "
          f"{float(DSQ)**.5:.5f} share a colour -- A FORCED PAIR."
          f"  [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump((T, [P[v] for v in T], DSQ),
                open(SC + f"BOTH_{M}_{DSQ}.pkl".replace("/", "_"), "wb"))
    break
print("DONE", flush=True)
