"""Circular ascent: grow the graph at the ratio where it is actually tight.

Every growth run in this project scored a candidate point by how many
5-colourings it kills, and every one of them stalled with a kill rate
decaying towards nothing.  The reason is now measurable.  De Grey's G has
chi_c well below 5 -- it maps to K(9/2) -- so asking about 5-colourings asks
a question with more than half a colour of slack in it.  Almost no single
point can kill a 5-colouring of a graph that loose, so the score carried no
information and the walk was blind.

The fix is to score at the ratio where the graph has no slack at all.  If
chi_c(H) = p/q then H maps to K(p/q) and to nothing below it: every
homomorphism to K(p/q) is on the boundary, and a point that kills one is
doing real work.  Push until no homomorphism to K(p/q) survives -- that is
an UNSAT proof that chi_c has risen -- then move to the next ratio up and
repeat.  The ladder of ratios between 4 and 5 is finite, its top is 5, and
passing 5 is exactly chi >= 6.

Unlike the integer chi, which answered "5" for every graph here and ranked
none of them, this objective moves.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
UNIV = sys.argv[1]                       # pkl holding the candidate pool
START = sys.argv[2]                      # "G" or a pkl for the seed
RATIO = Fr(*map(int, sys.argv[3].split("/")))
STEPS = int(sys.argv[4]) if len(sys.argv) > 4 else 400
SAMP = int(sys.argv[5]) if len(sys.argv) > 5 else 120
CAND = int(sys.argv[6]) if len(sys.argv) > 6 else 300
BATCH = int(sys.argv[7]) if len(sys.argv) > 7 else 25
t0 = time.time()

pool = pickle.load(open(SC + UNIV, "rb"))
seed = (build_G(K, as_graph=False) if START == "G"
        else pickle.load(open(SC + START, "rb")))
key = {}
allp = []
for P in (seed, pool):
    for x in P:
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
print(f"pool {len(allp)} points, seed {sum(inH)}, "
      f"{sum(len(v) for v in adj.values())//2} edges  "
      f"[{time.time()-t0:.0f}s]", flush=True)

p, q = RATIO.numerator, RATIO.denominator
W = [sum(1 << ((i + d) % p) for d in range(-(q - 1), q)) for i in range(p)]
FULL = (1 << p) - 1
print(f"ratio {RATIO} = {float(RATIO):.4f}; each neighbour forbids "
      f"{2*q-1} of {p} positions, so {-(-p//(2*q-1))} neighbours can block"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def clauses(verts):
    idx = {v: i for i, v in enumerate(verts)}
    cls = [[1 + i * p + j for j in range(p)] for i in range(len(verts))]
    for i in range(len(verts)):
        for j in range(p):
            for j2 in range(j + 1, p):
                cls.append([-(1 + i * p + j), -(1 + i * p + j2)])
    for v in verts:
        for u in adj[v]:
            if u in idx and u > v:
                for j in range(p):
                    for d in range(-(q - 1), q):
                        cls.append([-(1 + idx[v] * p + j),
                                    -(1 + idx[u] * p + (i2 := (j + d) % p))])
    return idx, cls


def homs(verts, want):
    """Sample distinct homomorphisms to K(p/q); [] means none exist."""
    idx, cls = clauses(verts)
    s = Solver(name="cd15", bootstrap_with=cls)
    out = []
    while len(out) < want and s.solve():
        m = s.get_model()
        phi = {}
        for v, i in idx.items():
            for j in range(p):
                if m[i * p + j] > 0:
                    phi[v] = j
                    break
        out.append(phi)
        blk = [-(1 + idx[v] * p + phi[v]) for v in random.sample(
            list(idx), min(len(idx), 40))]
        s.add_clause(blk)
    s.delete()
    return out


verts = [i for i in range(len(allp)) if inH[i]]
for step in range(STEPS):
    H = homs(verts, SAMP)
    if not H:
        print(f"\n*** no homomorphism to K({p}/{q}) at step {step}: "
              f"chi_c > {float(RATIO):.4f} on {len(verts)} points ***",
              flush=True)
        pickle.dump([allp[i] for i in verts],
                    open(SC + f"ascent_{p}_{q}.pkl", "wb"))
        break
    vs = set(verts)
    scored = []
    cands = [i for i in range(len(allp))
             if not inH[i] and len(adj[i] & vs) >= 3]
    for x in random.sample(cands, min(CAND, len(cands))):
        nb = list(adj[x] & vs)
        sc = 0
        for phi in H:
            m = 0
            for u in nb:
                m |= W[phi[u]]
                if m == FULL:
                    sc += 1
                    break
        scored.append((sc, x))
    if not scored:
        print("no candidate with three neighbours left", flush=True)
        break
    # A point that kills every sampled homomorphism is not rare at the tight
    # ratio -- at 9/2 three neighbours already suffice -- so adding one at a
    # time and resampling throws the sample away.  Take the whole top band.
    scored.sort(reverse=True)
    top = scored[0][0]
    batch = [x for sc, x in scored if sc == top][:BATCH]
    for x in batch:
        inH[x] = True
        verts.append(x)
    print(f"   step {step}: +{len(batch)} -> {len(verts)} points, "
          f"each killing {top}/{len(H)}, {len(cands)} candidates"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if step % 25 == 24:
        pickle.dump([allp[i] for i in verts],
                    open(SC + f"ascent_{p}_{q}.pkl", "wb"))
print("DONE", flush=True)
