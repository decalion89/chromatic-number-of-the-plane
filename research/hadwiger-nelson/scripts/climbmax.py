"""Hill-climb on one pair: which union raises ITS agreement?

G at five colours has 104 non-edge pairs agreeing in at least seventeen of
thirty-two sampled colourings, the best in twenty-one -- 0.66 against a chance
of 0.20, which is 3.30 times over.  Sa at four sat at 0.84 against 0.25, 3.36
times over, and one rotation took it to 1.00.

The rotation does not have to be about the pair's midpoint.  That is only the
symmetric choice, and it happens to be unavailable here: of G's top sixty
pairs, one sits on a rational closable ring about its own midpoint, and that
ring is D = 1, whose rotation lowers the agreement rather than raising it.

Any union that raises the pair's agreement will do.  So scan the pivots and
the closable rings, build G u rho(G) for each, and measure what happens to
THAT pair -- not to the graph in general.  Every scan in this file has been
undirected; this one has a target.
"""
import sys, time, random, pickle
from collections import Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SAMPLES, k = 32, 5
SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
P = build_G(F, as_graph=False)
g = build_graph(P)
n = g.n
E = set((min(a, b), max(a, b)) for a, b in g.edges())
GE = sorted(E)
ONE = F.rational(1)
zf0 = [(float(p.x), float(p.y)) for p in P]


def sample(pts, edges, seed=2718, samples=SAMPLES):
    m = len(pts)
    cls = [[1 + v * k + c for c in range(k)] for v in range(m)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None, None
    rng, cols = random.Random(seed), []
    for s in range(samples):
        sv.set_phases([-(1 + w) if rng.random() < .05 else (1 + w)
                       for w in range(m * k)])
        sv.solve()
        mo = sv.get_model()
        cols.append([next(c for c in range(k) if mo[w * k + c] > 0)
                     for w in range(m)])
    return np.array(cols, dtype=np.int8), sv


C, sv = sample(P, GE)
sv.delete()
targets = []
for i in range(n - 1):
    agree = (C[:, i + 1:] == C[:, i:i + 1]).sum(axis=0)
    for off in np.nonzero(agree >= 19)[0]:
        j = int(off) + i + 1
        if (i, j) not in E:
            targets.append((int(agree[off]), i, j))
targets.sort(reverse=True)
print(f"targets: {len(targets)} pairs at >= 19/{SAMPLES}; best "
      f"{[(a, i, j) for a, i, j, in targets[:5]]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
A0, I, J = targets[0]

# Every pivot with its closable rings, as in the undirected scans.
jobs = []
for v in range(n):
    px, py = zf0[v]
    rings = Counter()
    for q, (a, b) in enumerate(zf0):
        if q != v:
            rings[round((a - px) ** 2 + (b - py) ** 2, 9)] += 1
    for r2f, cnt in rings.items():
        D = Fr(round(r2f * 55440), 55440)
        if abs(float(D) - r2f) > 1e-8 or D < Fr(1, 4) or D == 1:
            continue
        if closable_distance(D):
            jobs.append((cnt, v, D))
jobs.sort(reverse=True)
print(f"{len(jobs)} (pivot, ring) candidates  [{time.time()-t0:.0f}s]",
      flush=True)

rng = random.Random(97)
rng.shuffle(jobs)
best = (A0, None)
for nth, (cnt, v, D) in enumerate(jobs[:400]):
    rho = rotation_joining(D, F).about(P[v])
    U, idx = list(P), {p: q for q, p in enumerate(P)}
    for p in P:
        q = rho(p)
        if q not in idx:
            idx[q] = len(U)
            U.append(q)
    zf = [(float(p.x), float(p.y)) for p in U]
    ed = set(GE)
    img = [idx[rho(p)] for p in P]
    for x, y in GE:
        u, w = img[x], img[y]
        ed.add((min(u, w), max(u, w)))
    cell = {}
    for q, (x, y) in enumerate(zf):
        cell.setdefault((int(x // 1), int(y // 1)), []).append(q)
    cross = 0
    for q in range(n, len(U)):
        x, y = zf[q]
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for w in cell.get((int(x // 1) + dx, int(y // 1) + dy), ()):
                    if w >= q or abs((x - zf[w][0]) ** 2
                                     + (y - zf[w][1]) ** 2 - 1) > 1e-7:
                        continue
                    if ((U[w].x - U[q].x) ** 2 + (U[w].y - U[q].y) ** 2
                            == ONE) and (w, q) not in ed:
                        ed.add((w, q))
                        cross += 1
    if cross == 0:
        continue
    CU, sv2 = sample(U, sorted(ed))
    if CU is None:
        print(f"  *** pivot {v} D={D}: NOT 5-COLOURABLE ***", flush=True)
        break
    # The maximum over ALL non-edge pairs, not one chosen in advance.  A
    # union that pushes some OTHER pair to the ceiling is exactly as good,
    # and following a fixed target would walk past it.  The reduction is one
    # numpy pass and costs nothing.
    m = len(U)
    es = set(ed)
    top, ti, tj = 0, -1, -1
    for x in range(m - 1):
        ag = (CU[:, x + 1:] == CU[:, x:x + 1]).sum(axis=0)
        y = int(ag.argmax())
        if int(ag[y]) > top and (x, y + x + 1) not in es:
            top, ti, tj = int(ag[y]), x, y + x + 1
    here = int((CU[:, I] == CU[:, J]).sum())
    if top > best[0]:
        forced = not sv2.solve(assumptions=[1 + ti * k, -(1 + tj * k)])
        best = (top, (v, str(D), cross, ti, tj, forced))
        print(f"  pivot {v} D={D}: best pair {ti},{tj} at {top}/{SAMPLES} "
              f"(target {I},{J} at {here}); {cross} cross"
              f"{'  *** FORCED ***' if forced else ''}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if forced:
            with open(SC + "climbpair_forced.pkl", "wb") as fh:
                pickle.dump((v, str(D), ti, tj,
                             [(str(p.x), str(p.y)) for p in U], sorted(ed)),
                            fh)
            sv2.delete()
            break
    sv2.delete()
    if nth % 50 == 0:
        print(f"  ... {nth}/400, best {best[0]}/{SAMPLES} at {best[1]}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"DONE: best {best[0]}/{SAMPLES} from {A0} at {best[1]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
