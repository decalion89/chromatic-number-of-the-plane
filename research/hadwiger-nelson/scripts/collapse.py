"""Collapse the two-branch disjunction to a single distance.

S is capped: every proper 5-colouring of G shows at most four colours on it,
so some colour class meets S in at least ceil(|S|/4) points, and that class is
independent.  No independent class of that size avoids both distance sqrt3 and
distance 2, which gives a disjunction over two distances.  A disjunction over
ONE distance would be a forced pair, and a forced pair at a closable distance
spindles into a sixth colour.

Any subset of a capped set is capped -- a colouring showing five colours on
the subset would show five on the whole -- so every T inside S inherits
palette at most 4 and carries its own pigeonhole, ceil(|T|/4).  That makes the
collapse a search inside S rather than a fresh search in G.

Shrinking T cuts the forced class as well as the escape routes, so the two
effects pull against each other and the question is which wins.  Handled the
standard way: find an escaping class, delete one of its points, repeat.  When
no escape survives, the single distance blocks and the disjunction is a forced
pair; when T runs down without that happening, it does not.

Run for each closable distance present, separately.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound, closable_distance
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
TARGET = int(sys.argv[2]) if len(sys.argv) > 2 else 4
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
Eset = set(E)
n = len(P)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
S, _ = pickle.load(open(SC + f"grown_{k}_{TARGET}.pkl", "rb"))
S = sorted(S)
print(f"S: {len(S)} points, palette "
      f"{ring_palette_bound(cls, n*k, S, k)} of {k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
dm, D2 = b.dim, b.D * b.D
pair_d = {}
for i, u in enumerate(S):
    for j in range(i + 1, len(S)):
        v = S[j]
        d = r[u] - r[v]
        sq = (b._field_square(d[None, :dm]) + b._field_square(d[None, dm:]))[0]
        if not any(sq[1:]):
            pair_d[(u, v)] = Fr(int(sq[0]), D2)


def escape(T, D):
    """An independent ceil(|T|/4)-subset of T with no pair at distance^2 D."""
    idx = {v: i for i, v in enumerate(T)}
    m = len(T)
    q = -(-m // TARGET)
    pool = IDPool(start_from=m + 1)
    f = list(CardEnc.atleast(lits=list(range(1, m + 1)), bound=q, vpool=pool,
                             encoding=EncType.seqcounter))
    for a in range(m):
        for c in range(a + 1, m):
            u, v = T[a], T[c]
            key = (min(u, v), max(u, v))
            if key in Eset or pair_d.get(key) == D:
                f.append([-(a + 1), -(c + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    out = [T[i] for i in range(m) if s.get_model()[i] > 0] if ok else None
    s.delete()
    return out, q


rng = random.Random(1049)
for D in (Fr(3), Fr(4)):
    print(f"\n--- trying to make distance^2 {D} block alone "
          f"(d = {float(D)**.5:.5f}, closable "
          f"{bool(closable_distance(D))})", flush=True)
    T = list(S)
    for step in range(400):
        esc, q = escape(T, D)
        if esc is None:
            pal = ring_palette_bound(cls, n * k, T, k)
            print(f"*** BLOCKS ALONE: {len(T)} points, palette {pal}, forced "
                  f"class >= {q}, and no such class avoids distance^2 {D}.  "
                  f"So in EVERY 5-colouring of G some pair at distance "
                  f"{float(D)**.5:.5f} is monochromatic -- A FORCED PAIR."
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            with open(SC + f"forcedpair_{D}.pkl", "wb") as fh:
                pickle.dump((T, [P[v] for v in T], D), fh)
            break
        # drop a point of the escaping class, preferring one in few D-pairs
        def dcount(v):
            return sum(1 for w in T if w != v
                       and pair_d.get((min(v, w), max(v, w))) == D)
        victim = min(esc, key=lambda v: (dcount(v), rng.random()))
        T = [v for v in T if v != victim]
        if len(T) < 8:
            print(f"   ran down to {len(T)} points without blocking"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            break
        if step % 10 == 0:
            print(f"   step {step}: {len(T)} points, forced class {q}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
