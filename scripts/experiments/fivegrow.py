"""Grow by five-cycles, which is what the orientation analysis asks for.

Every growth run before this scored a point by how many sampled colourings
or homomorphisms it killed, and the sample is a vanishing part of the space,
so the score carried no information.  The Minty reading names a structural
target instead: an orientation survives at 9/2 only if EVERY 5-cycle gets
two edges on its minority side, a lone 5-cycle always can, and overlapping
ones sharing edges cannot always.  So the thing to add is 5-cycles that
overlap what is already there.

The score needs no solver at all.  A 5-cycle through a new point x runs
x - a - b - c - d - x with a and d neighbours of x, so counting them is
counting 3-edge paths between pairs of x's neighbours.  That is cheap,
deterministic, and measures the structure rather than a draw from it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/hn/")
POOL, SEED = sys.argv[1], sys.argv[2]
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 400
CAND = int(sys.argv[4]) if len(sys.argv) > 4 else 400
BATCH = int(sys.argv[5]) if len(sys.argv) > 5 else 20
OUT = sys.argv[6] if len(sys.argv) > 6 else "fivegrow.pkl"
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
print(f"pool {len(allp)}, seed {len(verts)}  [{time.time()-t0:.0f}s]",
      flush=True)


def five_through(x, vs):
    """5-cycles x-a-b-c-d-x with a,b,c,d all inside vs."""
    nb = sorted(adj[x] & vs)
    if len(nb) < 2:
        return 0
    nbs = set(nb)
    tot = 0
    for a in nb:
        for bb in adj[a] & vs:
            if bb == x:
                continue
            for c in adj[bb] & vs:
                if c in (x, a):
                    continue
                tot += len((adj[c] & nbs) - {a, bb, x})
    return tot // 2


def density(vs):
    m = sum(1 for v in vs for u in adj[v] if u in vs) // 2
    tot = sum(five_through(v, vs) for v in random.sample(sorted(vs), 60))
    return m, tot / 60


random.seed(9)
m0, d0 = density(verts)
print(f"seed: {m0} edges, {d0:.1f} five-cycles per vertex"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for step in range(STEPS):
    cands = [i for i in range(len(allp))
             if not inH[i] and len(adj[i] & verts) >= 2]
    if not cands:
        break
    scored = []
    for x in random.sample(cands, min(CAND, len(cands))):
        scored.append((five_through(x, verts), x))
    scored.sort(reverse=True)
    batch = [x for sc, x in scored[:BATCH] if sc > 0]
    if not batch:
        print("no candidate creates a five-cycle", flush=True)
        break
    for x in batch:
        inH[x] = True
        verts.add(x)
    if step % 5 == 0 or step < 3:
        m = sum(1 for v in verts for u in adj[v] if u in verts) // 2
        print(f"   step {step}: +{len(batch)} -> {len(verts)} pts, "
              f"{m} edges, best adds {scored[0][0]} five-cycles"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if step % 10 == 9:
        pickle.dump([allp[i] for i in sorted(verts)], open(SC + OUT, "wb"))
pickle.dump([allp[i] for i in sorted(verts)], open(SC + OUT, "wb"))
m1, d1 = density(verts)
print(f"\nfinal: {len(verts)} points, {m1} edges, "
      f"{d1:.1f} five-cycles per vertex (was {d0:.1f})"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
