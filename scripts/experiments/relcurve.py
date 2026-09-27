"""Does the colour relation switch on with size?  Put points in the gap.

At four colours Sa is fifty-seven times the minimum for chi = 4 and its
relation sits twenty times above chance.  At five colours G is three times the
minimum and sits at chance.  Between three and fifty-seven, nothing has been
measured at five colours, and a single point either side is not a curve.

So measure the ratio at five colours for a series of graphs of growing size,
at a FIXED sample count so the numbers compare: G itself, a pivot union of G,
a translate union, and two depths of the translate stack.  If the ratio climbs
with size the hypothesis lives and says how far it has to go; if it stays at
one, size is not the variable and the whole family is out of reach.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import agreeing_pairs
from pysat.solvers import Solver

SC = ("/tmp/hn/")
SAMPLES, k, MIN5 = 40, 5, 500
t0 = time.time()


def ratio(name, pts, verify=1500):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT 5-colourable", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(4096), []
    for s in range(SAMPLES):
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
        msk = (C[:, i + 1:] != C[:, i:i + 1]).all(axis=0)
        for off in np.nonzero(msk)[0]:
            j = int(off) + i + 1
            if (i, j) not in E:
                diff.append((i, j))
    pr = n * (n - 1) // 2
    ch = pr * ((k - 1) / k) ** SAMPLES
    fs = [(i, j) for i, j in same
          if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    fd = [(i, j) for i, j in diff[:verify]
          if not sv.solve(assumptions=[1 + i * k, 1 + j * k])]
    sv.delete()
    print(f"  {name}: {n} points ({n/MIN5:.1f}x the minimum), {pr} pairs, "
          f"chance {ch:.1f}, {len(diff)} candidates, ratio "
          f"{len(diff)/max(ch,1e-9):.2f}; VERIFIED {len(fs)} same, {len(fd)} "
          f"different of {min(len(diff), verify)} checked  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return len(diff) / max(ch, 1e-9)


P = build_G(F, as_graph=False)
ratio("G", P)

# A pivot union: the best-crossing rotation the scan found, D = 4 at pivot 14.
g0 = build_graph(P)
rho = rotation_joining(Fr(4), F).about(list(g0.vertices)[14])
U = list(P)
seen = set(P)
for p in P:
    q = rho(p)
    if q not in seen:
        seen.add(q)
        U.append(q)
ratio("G u rho_4(G) about vertex 14", U)

# Translate stack, depths 1 and 2, from the census leaders.
with open(SC + "tcensus.pkl", "rb") as fh:
    CEN = pickle.load(fh)
RAD = []
for mask in range(16):
    pr = 1
    for b in range(4):
        if mask >> b & 1:
            pr *= (3, 5, 7, 11)[b]
    RAD.append(pr ** 0.5)
cur = list(P)
curset = set(P)
for depth in (1, 2):
    t = CEN[depth - 1][1]
    tx = F.zero()
    ty = F.zero()
    for m in range(16):
        if t[0][m]:
            tx = tx + F.sqrt(int(round(RAD[m] ** 2))) * F.rational(t[0][m]) \
                if RAD[m] != 1 else tx + F.rational(t[0][m])
        if t[1][m]:
            ty = ty + F.sqrt(int(round(RAD[m] ** 2))) * F.rational(t[1][m]) \
                if RAD[m] != 1 else ty + F.rational(t[1][m])
    nxt = list(cur)
    for p in cur:
        q = Point(p.x + tx, p.y + ty)
        if q not in curset:
            curset.add(q)
            nxt.append(q)
    cur = nxt
    ratio(f"translate stack depth {depth}", cur)
