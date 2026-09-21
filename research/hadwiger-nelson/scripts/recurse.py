"""Does spindling shrink the obstruction?  Iterate and count.

Spindling the witness over all 161 forbidden pairs gives 5363 points and it
5-colours.  That confirms by construction what was only objected to in the
abstract: rotated copies do not consume a disjunction, because each copy is
free to have a DIFFERENT pair of the list monochromatic.

But the union is a new graph, and it has its own obstruction.  The question
that matters is whether spindling makes the obstruction SMALLER.  If each
round cuts the number of pairs, iterating drives it to one -- and a single
forced pair is de Grey's hypothesis exactly, which the ordinary two-copy
spindle does consume.  If the count holds or grows, the recursion is a
treadmill and that is worth knowing in one number rather than guessed at.

So per round: build the union, find its minimal obstruction the same way --
cumulative classes, then greedy removal of pairs and vertices -- and report
the count.  The sequence is the answer.
"""
import sys, time, json, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.field import Field
from hn.geometry import Point, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
doc = json.load(open("/home/user/darwin-50/research/hadwiger-nelson/"
                     "data/witness_five.json"))
K = Field(tuple(doc["field"]["generators"]))
P = [Point(K.element([Fr(c) for c in x]), K.element([Fr(c) for c in y]))
     for x, y in doc["points"]]
FORB = [(Fr(D), (a, b)) for a, b, D in doc["forbidden_pairs"]]


def graph_of(pts):
    basis = IntBasis.covering(pts)
    rows = basis.rows(pts)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    return basis, rows, E


def classes_of(pts, basis, rows, E):
    dim, D2 = basis.dim, basis.D * basis.D
    Eset = set(E)
    byd = defaultdict(list)
    for i in range(len(pts) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in Eset:
                v = int(sq[off, 0])
                if v:
                    dd = Fr(v, D2)
                    if closable_distance(dd):
                        byd[dd].append((i, j))
    return byd


def colours(n, E, prs):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    for _, (a, b) in prs:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


def obstruction(pts, cap=400):
    basis, rows, E = graph_of(pts)
    n = len(pts)
    byd = classes_of(pts, basis, rows, E)
    clo = sorted(byd, key=lambda d: len(byd[d]))
    print(f"      {len(clo)} closable classes, sizes "
          f"{[len(byd[d]) for d in clo[:8]]}...{[len(byd[d]) for d in clo[-4:]]}",
          flush=True)
    cum = []
    for D in clo[:cap]:
        cum += [(D, p) for p in byd[D]]
        if not colours(n, E, cum):
            break
    else:
        return None, n, len(E), None
    # greedy pair removal, a few orders
    best = list(cum)
    for t in range(3):
        cur = list(cum)
        rng = random.Random(90 + t)
        while True:
            snap = (sorted(cur, key=lambda x: -float(x[0])) if t == 0
                    else rng.sample(cur, len(cur)))
            dropped = 0
            for target in snap:
                if target not in cur:
                    continue
                trial = [x for x in cur if x is not target]
                if colours(n, E, trial):
                    continue
                cur = trial
                dropped += 1
            if dropped == 0:
                break
        if len(cur) < len(best):
            best = cur
    return best, n, len(E), len(set(D for D, _ in best))


cur_pts, cur_forb = P, FORB
for rnd in range(1, 5):
    c = Counter(str(D) for D, _ in cur_forb)
    print(f"\nround {rnd}: {len(cur_pts)} points, {len(cur_forb)} pairs over "
          f"{len(c)} classes {dict(sorted(c.items()))}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    seen, U = set(cur_pts), list(cur_pts)
    for D, (a, b) in cur_forb:
        r = rotation_joining(D, K).about(cur_pts[a])
        for p in cur_pts:
            q = r(p)
            if q not in seen:
                seen.add(q)
                U.append(q)
    basis, rows, E = graph_of(U)
    n = len(U)
    ok = colours(n, E, [])
    print(f"   spindled: {n} points, {len(E)} edges -> "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
    if n > 30000:
        print("   too big to continue", flush=True)
        break
    nf, nn, ne, ncl = obstruction(U)
    if nf is None:
        print(f"   no obstruction found within the class cap", flush=True)
        break
    print(f"   its own obstruction: {len(nf)} pairs over {ncl} classes "
          f"(was {len(cur_forb)} over {len(c)})  [{time.time()-t0:.0f}s]",
          flush=True)
    cur_pts, cur_forb = U, nf
print("DONE", flush=True)
