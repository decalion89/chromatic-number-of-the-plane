"""Filter by sampling, decide by SAT -- and the whole census becomes reachable.

Separating all 44011 candidate pairs of a union costs the solver ~700000
conflicts, twenty minutes a translate.  At 399 translates that is a week.  But
a pair is forced only if it agrees in EVERY proper 5-colouring, so a single
colouring in which it differs refutes it outright.  Sample colourings instead
and the candidate set collapses geometrically: a random colouring kills four
fifths of the pairs, the next four fifths of those, and after seven samples
essentially nothing is left.  Only the survivors get a SAT check, and only
those checks can be expensive.

The filter is sound in the direction that matters.  A pair that differs in a
sampled colouring is definitively not forced -- the colouring is the witness.
So the sampling never hides a forced pair; it only refuses to confirm one, and
the SAT check does the confirming.

Diversity comes from assuming a handful of random vertices into random
colours.  Phase saving would otherwise hand back the same colouring every
time, and identical samples filter nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from msqrt import madd, msub, mmul
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

GENS = (3, 5, 7, 11)
DIM = 16
SAMPLES = 120
STALL = 25
SC = ("/tmp/hn/")
random.seed(20260920)
t0 = time.time()
pts = build_G(F, as_graph=False)
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
zf = [(float(p.x), float(p.y)) for p in pts]
GE = list(build_graph(pts).edges())
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)
RAD = []
for m in range(DIM):
    pr = 1
    for b in range(4):
        if m >> b & 1:
            pr *= GENS[b]
    RAD.append(pr ** 0.5)


def fl(c):
    return sum(float(x) * r for x, r in zip(c, RAD))


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


with open(SC + "tcensus.pkl", "rb") as fh:
    CEN = pickle.load(fh)
print(f"census: {len(CEN)} translates, leaders "
      f"{[c for c, _ in CEN[:8]]}  [{time.time()-t0:.0f}s]", flush=True)

Pidx = {p: i for i, p in enumerate(P)}
best_surv = -1
for ti, (ccount, t) in enumerate(CEN):
    tz = (fl(t[0]), fl(t[1]))
    allp, allz = list(P), list(zf)
    idx = dict(Pidx)
    second = []
    for p, z in zip(P, zf):
        q = (madd(p[0], t[0]), madd(p[1], t[1]))
        if q in idx:
            second.append(idx[q])
        else:
            idx[q] = len(allp)
            second.append(len(allp))
            allp.append(q)
            allz.append((z[0] + tz[0], z[1] + tz[1]))
    n0, n = len(P), len(allp)
    cell = {}
    for i in range(n0):
        a, b = allz[i]
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    E = set()
    for a0, b0 in GE:
        E.add((min(a0, b0), max(a0, b0)))
        x, y = second[a0], second[b0]
        E.add((min(x, y), max(x, y)))
    cross = 0
    for j in range(n0, n):
        a, b = allz[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if abs((a - allz[i][0]) ** 2
                           + (b - allz[i][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = (msub(allp[i][0], allp[j][0]),
                         msub(allp[i][1], allp[j][1]))
                    if norm2(d) == ONE and (i, j) not in E:
                        E.add((i, j))
                        cross += 1
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in sorted(E):
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** translate {ti}: {n} pts, {cross} cross -- NOT "
              f"5-COLOURABLE ***  [{time.time()-t0:.0f}s]", flush=True)
        with open(SC + "tfast_six.pkl", "wb") as fh:
            pickle.dump((ti, t, allp, sorted(E)), fh)
        sys.exit(0)
    cand = []
    for i in range(n):
        ai, bi = allz[i]
        for j in range(i + 1, n):
            v = (ai - allz[j][0]) ** 2 + (bi - allz[j][1]) ** 2
            if v > 36.0:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
                continue
            cand.append((i, j))
    surv = cand
    since, used = 0, 0
    for s in range(SAMPLES):
        for _try in range(6):
            vs = random.sample(range(n), 8)
            asm = [1 + v * 5 + random.randrange(5) for v in vs]
            if sv.solve(assumptions=asm):
                break
        else:
            sv.solve()
        m = sv.get_model()
        col = [0] * n
        for w in range(n):
            for c in range(5):
                if m[w * 5 + c] > 0:
                    col[w] = c
                    break
        used = s + 1
        before = len(surv)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        since = 0 if len(surv) < before else since + 1
        if not surv or since >= STALL:
            break
    hits = [(i, j) for i, j in surv
            if not sv.solve(assumptions=[1 + i * 5, -(1 + j * 5)])]
    sv.delete()
    tag = "" if len(surv) <= best_surv else "  <-- best"
    best_surv = max(best_surv, len(surv))
    print(f"  {ti}: {n} pts, {cross} cross, {len(cand)} pairs -> "
          f"{len(surv)} survive {used} samples, {len(hits)} FORCED"
          f"  [{time.time()-t0:.0f}s]{tag}", flush=True)
    if hits:
        with open(SC + "tfast_forced.pkl", "wb") as fh:
            pickle.dump((ti, t, allp, sorted(E), hits), fh)
        print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
              flush=True)
        break
