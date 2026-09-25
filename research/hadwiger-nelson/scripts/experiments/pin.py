"""Grow against the free fraction, the one objective with a calibrated target.

Every growth objective in this project was uncalibrated: kill counts with no
known target, saturation scores with no reference, five-cycle densities with
no threshold.  The free fraction has all three.  Take a proper colouring and
count the vertices whose neighbours leave them a spare colour; Sa at four
reaches ZERO of 397, which is where de Grey's construction works, and G at
five sits at 22.1 per cent.  The target is zero and something achieves it
one level down.

union3 reaches 18.9 per cent with twenty-five times G's points, but union3
was never GROWN for this -- it is a generic closure.  Here the score is the
objective itself: add the points that touch the most currently-free
vertices, since a free vertex loses its freedom only when a new neighbour
takes its spare colour.

No solver in the scoring loop.  One colouring per round, a count, and a
choice.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
POOL = sys.argv[1] if len(sys.argv) > 1 else "union3.pkl"
SEED = sys.argv[2] if len(sys.argv) > 2 else "G"
KC = int(sys.argv[3]) if len(sys.argv) > 3 else 5
STEPS = int(sys.argv[4]) if len(sys.argv) > 4 else 300
BATCH = int(sys.argv[5]) if len(sys.argv) > 5 else 30
t0 = time.time()
pool = pickle.load(open(SC + POOL, "rb"))
seed = (build_G(K, as_graph=False) if SEED == "G"
        else pickle.load(open(SC + SEED, "rb")))
key, allp = {}, []
for Q in (seed, pool):
    for x in Q:
        k = (tuple(x.x.c), tuple(x.y.c))
        if k not in key:
            key[k] = len(allp)
            allp.append(x)
inH = [False] * len(allp)
for x in seed:
    inH[key[(tuple(x.x.c), tuple(x.y.c))]] = True
b = IntBasis.covering(allp)
rows = b.rows(allp)
assert b.overflow_headroom(rows) < 1.0
adj = defaultdict(set)
for a, c in fast_edges_complete(b, rows):
    adj[a].add(c)
    adj[c].add(a)
verts = set(i for i in range(len(allp)) if inH[i])
print(f"pool {len(allp)}, seed {len(verts)}, k = {KC}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def colour_and_free(vs):
    """One proper colouring, and the set of vertices with a spare colour."""
    ks = sorted(vs)
    idx = {v: i for i, v in enumerate(ks)}
    cls = [[1 + i * KC + c for c in range(KC)] for i in range(len(ks))]
    for v in ks:
        for u in adj[v]:
            if u in idx and u > v:
                for c in range(KC):
                    cls.append([-(1 + idx[v] * KC + c),
                                -(1 + idx[u] * KC + c)])
    s = Solver(name="m22", bootstrap_with=cls)
    if not s.solve():
        s.delete()
        return None, None
    m = s.get_model()
    s.delete()
    col = {}
    for v in ks:
        for c in range(KC):
            if m[idx[v] * KC + c] > 0:
                col[v] = c
                break
    free = set()
    for v in ks:
        used = {col[u] for u in adj[v] if u in vs}
        if len(used | {col[v]}) < KC:
            free.add(v)
    return col, free


col, free = colour_and_free(verts)
print(f"seed: {len(free)} of {len(verts)} free = "
      f"{100*len(free)/len(verts):.1f}%  [{time.time()-t0:.0f}s]", flush=True)
random.seed(23)
best = len(free) / len(verts)
for step in range(STEPS):
    # Adding a point adds a vertex that may itself be free, so the fraction
    # can RISE even while each addition constrains others -- which is what
    # happened: 13.38 per cent at 7611 points, back to 16.14 at 10311.  A
    # candidate must be able to end up pinned itself, which needs at least k
    # neighbours already present.
    cands = [i for i in range(len(allp))
             if not inH[i] and (adj[i] & free)
             and len(adj[i] & verts) >= KC]
    if not cands:
        print("no candidate touches a free vertex", flush=True)
        break
    scored = sorted(((len(adj[i] & free), i) for i in
                     random.sample(cands, min(3000, len(cands)))),
                    reverse=True)
    batch = [i for sc, i in scored[:BATCH] if sc > 0]
    if not batch:
        break
    for x in batch:
        inH[x] = True
        verts.add(x)
    col, free = colour_and_free(verts)
    if col is None:
        print(f"   *** step {step}: {len(verts)} points and NOT "
              f"{KC}-COLOURABLE ***", flush=True)
        pickle.dump([allp[i] for i in sorted(verts)],
                    open(SC + "pinned.pkl", "wb"))
        break
    r = len(free) / len(verts)
    if r < best:
        best = r
        pickle.dump([allp[i] for i in sorted(verts)],
                    open(SC + "pinned.pkl", "wb"))
    if step % 5 == 0:
        print(f"   step {step}: +{len(batch)} -> {len(verts)} points, "
              f"{len(free)} free = {100*r:.2f}%  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\nbest free fraction {100*best:.2f}% (G starts at 22.1, union3 "
      f"reaches 18.9, Sa at four colours reaches 0.0)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
