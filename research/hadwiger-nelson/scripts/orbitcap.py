"""A capped set that is a union of orbits, so that its disjunction is symmetric.

The capped set found in G gives a real disjunction -- some pair at sqrt3 or at
2 is monochromatic in every 5-colouring -- and it cannot be used, for a reason
that took a while to see clearly.

A FORCED pair closes immediately: rotate a copy of the graph about u until v
and its image are a unit apart, and all three share a colour while two of them
are adjacent.  A DISJUNCTION does not, because in the rotated copy the
monochromatic pair may be a different one, and no contradiction follows.

de Grey's disjunction survives that because his three antipodal pairs are an
ORBIT of the hexagon's symmetry.  The rotation maps the disjunction onto
itself, so whichever pair is monochromatic, its image is again one of the
three, and the argument closes.  The bite is not a trick for collapsing a
disjunction in general -- it works because the disjunction is symmetric.

So the capped set has to be symmetric too.  Gc, the closure of G under the
order-twelve group about the centre with the most unit-distance neighbours, is
the carrier where that is possible: its points fall into group orbits, and any
union of orbits is mapped to itself by the rotation.  This asks the decision
procedure for a capped union of orbits, with the orbit as the unit of choice
rather than the point -- which also shrinks the search from eleven thousand
variables to a few hundred.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NORB = int(sys.argv[2]) if len(sys.argv) > 2 else 3
TARGET = int(sys.argv[3]) if len(sys.argv) > 3 else 4
SEEDN = int(sys.argv[4]) if len(sys.argv) > 4 else 300
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = pickle.load(open(SC + "wide0.pkl", "rb"))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
deg = defaultdict(int)
for a, c in E:
    deg[a] += 1
    deg[c] += 1
# The centre is the one the closure was built about -- G's anchor point --
# not whichever vertex happens to have the highest degree.  Getting this wrong
# makes the computed orbits meaningless, and it shows: under a group of order
# twelve every orbit has size dividing twelve, so sizes 2, 4 and 8 are a
# symptom, not data.
from hn.degrey import build_G
C = build_G(K, as_graph=False)[0]
ci = next(i for i, q in enumerate(P) if q == C)
print(f"Gc: {n} points, {len(E)} edges; closure centre is vertex {ci} at "
      f"({float(C.x):.4f},{float(C.y):.4f}), degree {deg[ci]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
rot = _rot60(K).about(C)
idx = {p: i for i, p in enumerate(P)}
orbit = [None] * n
orbs = []
for i, p in enumerate(P):
    if orbit[i] is not None:
        continue
    mem = set()
    for base in (p, Point(p.x, C.y + (C.y - p.y))):
        q = base
        for _ in range(6):
            j = idx.get(q)
            if j is not None:
                mem.add(j)
            q = rot(q)
    o = len(orbs)
    fresh = sorted(j for j in mem if orbit[j] is None)
    for j in fresh:
        orbit[j] = o
    orbs.append(sorted(mem))
    assert len(mem) in (1, 2, 3, 4, 6, 12), (len(mem), "not a group orbit")
sizes = defaultdict(int)
for o in orbs:
    sizes[len(o)] += 1
print(f"{len(orbs)} orbits, sizes {dict(sorted(sizes.items()))}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
big = [i for i, o in enumerate(orbs) if len(o) >= 6]
print(f"{len(big)} orbits of at least 6 points  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
col_solver = Solver(name="cd15", bootstrap_with=cls)
assert col_solver.solve(), "Gc is not 5-colourable?"
rng = random.Random(32452843)


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
print(f"seeded with {len(samples)} colourings; choosing {NORB} orbits"
      f"  [{time.time()-t0:.0f}s]", flush=True)
NV = len(big)
for rnd in range(20000):
    pool = IDPool(start_from=NV + 1)
    y = list(range(1, NV + 1))
    f = list(CardEnc.equals(lits=y, bound=NORB, vpool=pool,
                            encoding=EncType.seqcounter))
    for si, c in enumerate(samples):
        # a colour is missing from the chosen union iff it is missing from
        # every chosen orbit
        miss = [pool.id(("z", si, q)) for q in range(k)]
        f += list(CardEnc.atleast(lits=miss, bound=k - TARGET, vpool=pool,
                                  encoding=EncType.seqcounter))
        for q in range(k):
            for t, oi in enumerate(big):
                if any(c[v] == q for v in orbs[oi]):
                    f.append([-miss[q], -(t + 1)])
    s = Solver(name="cd15", bootstrap_with=f)
    ok = s.solve()
    mod = s.get_model() if ok else None
    s.delete()
    if not ok:
        print(f"\nUNSATISFIABLE after {len(samples)} colourings: no union of "
              f"{NORB} orbits of Gc is capped at {TARGET}.  An absence proof."
              f"  [{time.time()-t0:.0f}s]", flush=True)
        break
    chosen = [big[t] for t in range(NV) if mod[t] > 0]
    T = sorted(v for oi in chosen for v in orbs[oi])
    w = witness(T)
    if w is None:
        pal = ring_palette_bound(cls, n * k, T, k)
        print(f"\n*** SYMMETRIC CAPPED SET: {len(chosen)} orbits, {len(T)} "
              f"points, palette {pal} of {k} -- and it is mapped to itself by "
              f"the sixty-degree rotation about the centre"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        pickle.dump((T, [P[v] for v in T], chosen, C),
                    open(SC + f"orbitcap_{NORB}.pkl", "wb"))
        break
    samples.append(w)
    if rnd % 100 == 0 and rnd:
        print(f"   round {rnd}: {len(samples)} colourings"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
