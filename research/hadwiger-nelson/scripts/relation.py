"""Is G's colour relation at five colours trivial?  Both halves of it.

A graph constrains its colourings through two relations: pairs forced to
agree, and pairs forced to differ.  Edges give the second for free, so the
question is whether anything beyond them survives.

The agreement half is already answered -- 0 pairs across 123 unions and every
graph measured here.  This is the other half, and it is filtered the same way
but inverted: a pair that AGREES in some sampled colouring cannot be forced
apart, so the candidates are the pairs that differ in all of them.  At
twenty-four samples a free pair survives with probability (4/5)^24, about one
in two hundred, which leaves a few thousand out of a million -- and each of
those is one fast SAT call at five colours.

If both halves come back empty, G's colour relation is exactly its edge set:
no pair of non-adjacent points is constrained in either direction, in any
proper 5-colouring, and there is nothing for a rotation to work with.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()


def relation(name, pts, k, samples=24, seed=6180, verify=True):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(seed), []
    for s in range(samples):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    same = [(i, j) for i, j in agreeing_pairs(cols, colours=k)
            if (i, j) not in E]
    C = np.array(cols, dtype=np.int8)
    diff = []
    for i in range(n - 1):
        m = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
        for off in np.nonzero(m)[0]:
            j = int(off) + i + 1
            if (i, j) not in E:
                diff.append((i, j))
    print(f"  {name}: {n} points, {len(E)} edges, {samples} samples -> "
          f"{len(same)} candidate forced-same, {len(diff)} candidate "
          f"forced-different  [{time.time()-t0:.0f}s]", flush=True)
    if not verify:
        # At four colours a single forced-different proof runs to minutes and
        # there are fifteen hundred candidates.  The counts are the signal:
        # twenty-four samples leave a free pair with probability (4/5)^24,
        # about one in two hundred, so anything far above that is structure.
        # The chance rate is ((k-1)/k)^samples, not 0.8^samples: two random
        # vertices differ with probability 3/4 at four colours and 4/5 at
        # five.  Using 0.8 for both understates the four-colour contrast by a
        # factor of five.
        print(f"    (four colours: counts only; chance would leave about "
              f"{int(n*(n-1)/2*((k-1)/k)**samples)})", flush=True)
        sv.delete()
        return
    fs = [(i, j) for i, j in same
          if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    fd = [(i, j) for i, j in diff
          if not sv.solve(assumptions=[1 + i * k, 1 + j * k])]
    sv.delete()
    print(f"    VERIFIED: {len(fs)} forced-same, {len(fd)} forced-different "
          f"beyond the edges  [{time.time()-t0:.0f}s]", flush=True)
    if fs or fd:
        with open(SC + f"relation_{name.replace(' ', '_')}.pkl", "wb") as fh:
            pickle.dump((name, k, fs, fd), fh)


relation("G at five", build_G(F, as_graph=False), 5)
relation("Y at five", build_Y(F), 5)
relation("Sa at five", build_Sa(F), 5)
relation("Sa at four (control)", build_Sa(F), 4, verify=False)
relation("Y at four (control)", build_Y(F), 4, verify=False)
