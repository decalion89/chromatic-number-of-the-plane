"""Spindle the witness and just ask, since it is small enough to build.

The witness says: in every proper 5-colouring, one of 161 pairs is
monochromatic.  Turning that into chi >= 6 by argument needs the disjunction
consumed, and the two-copy spindle does not consume a disjunction -- the
rotated copy may have a DIFFERENT pair of the list monochromatic, and no
contradiction closes.  That objection stands.

It does not stop the construction.  For each forbidden pair (a, b) at squared
distance D, the spindle rotation for D about the end a sends b to a point one
away from it; union the witness with all of those rotated copies and ask the
solver outright.  If the union fails to colour, that IS chi >= 6, with no
argument about disjunctions needed -- the SAT answer is the proof, and the
object is explicit.

On G this was unthinkable: 161 copies of 1581 points is a quarter of a
million.  On the 63-point witness it is ten thousand, which is one call.

Tried in stages -- one class of pairs at a time, then all of them -- so that
a smaller union gets its chance before the big one, and so the sizes and
verdicts are recorded either way.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.field import Field
from hn.geometry import Point, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 5
t0 = time.time()
doc = json.load(open("/home/user/darwin-50/research/hadwiger-nelson/"
                     "data/witness_five.json"))
K = Field(tuple(doc["field"]["generators"]))
P0 = [Point(K.element([Fr(c) for c in x]), K.element([Fr(c) for c in y]))
      for x, y in doc["points"]]
FORB = [(Fr(D), (a, b)) for a, b, D in doc["forbidden_pairs"]]
byd = defaultdict(list)
for D, p in FORB:
    byd[D].append(p)
print(f"witness: {len(P0)} points, {len(FORB)} forbidden pairs over "
      f"{len(byd)} classes  [{time.time()-t0:.0f}s]", flush=True)


def spindle_rotation(D):
    """The rotation taking a point at squared distance D from the centre to
    one exactly 1 away from where it was."""
    return rotation_joining(D, K)


def build(classes):
    seen, U = set(P0), list(P0)
    copies = 0
    for D in classes:
        rot = spindle_rotation(D)
        for a, b in byd[D]:
            r = rot.about(P0[a])
            copies += 1
            for p in P0:
                q = r(p)
                if q not in seen:
                    seen.add(q)
                    U.append(q)
    return U, copies


def solve(U):
    basis = IntBasis.covering(U)
    rows = basis.rows(U)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    n = len(U)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return n, len(E), ok


# sanity: the spindle rotation really moves the far end by one
for D in sorted(byd):
    a, b = byd[D][0]
    r = spindle_rotation(D).about(P0[a])
    moved = P0[b].dist2(r(P0[b]))
    print(f"   D={D}: the far end moves by {moved} (want 1)", flush=True)

for classes in ([[D] for D in sorted(byd)] + [sorted(byd)]):
    U, copies = build(classes)
    if len(U) > 60000:
        print(f"  classes {[str(d) for d in classes]}: {len(U)} points, "
              f"too big  [{time.time()-t0:.0f}s]", flush=True)
        continue
    n, m, ok = solve(U)
    tag = ",".join(str(d) for d in classes)
    print(f"  classes [{tag}]: {copies} copies -> {n} points, {m} edges -> "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
