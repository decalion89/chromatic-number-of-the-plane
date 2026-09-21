"""Buy contact instead of structure, and ask the cap again.

De Grey's bite joins Sa to its image by SIX edges.  Sweeping every rational
ring the field can join shows that was a structural choice, not a maximal one:
D = 5/9 and D = 5/11 give 108 crossing edges, eighteen times as many, and
every rotation about the origin commutes with the sixty-degree rotation so
every union stays dihedrally symmetric with its ring structure intact.

So take them all.  Sa together with its images under the highest-contact
rotations in both directions is the best-connected symmetric object this
family offers, at about the scale the measured factor of ten predicts for five
colours.  Then ask the questions that matter, in increasing strength:

  - does it colour with five at all;
  - is any ring about the origin CAPPED below five, which is what de Grey's
    construction really rests on;
  - does any ring carry the gateway, its antipodal corollary.

Built cumulatively so the trend is visible rather than one number at the end.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, doubly_usable_ring
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
Sa = build_Sa(K)
ZERO = Point(K.zero(), K.zero())
ORDER = [Fr(5, 9), Fr(5, 11), Fr(15, 16), Fr(3, 5), Fr(3, 7), Fr(15, 11),
         Fr(4), Fr(25)]
rots = []
for D in ORDER:
    rot = rotation_joining(D, K)
    rots.append((D, "+", rot))
    rots.append((D, "-", Rotation(rot.cos, -rot.sin)))
print(f"Sa {len(Sa)} points; {len(rots)} rotations, k={k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def study(U, label):
    b = IntBasis.covering(U)
    r = b.rows(U)
    head = b.overflow_headroom(r)
    if head >= 1.0:
        print(f"  {label}: headroom {head:.2f} -- NOT CERTIFIED, stop",
              flush=True)
        return None
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(U)
    zi = U.index(ZERO)
    d = r - r[zi]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), d2)
        if v:
            grp[v].append(int(off))
    rings = [(D, mem) for D, mem in grp.items()
             if len(mem) >= 4 and closable_distance(D)]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    nv = n * k
    s0 = Solver(name="cd15", bootstrap_with=cls)
    base = s0.solve()
    s0.delete()
    if not base:
        print(f"  *** {label}: {n} pts NOT {k}-COLOURABLE -- chi >= {k+1} ***",
              flush=True)
        return "WIN"
    # the cap screen: one call per ring
    ind = [nv + 1 + c for c in range(k)]
    enc = CardEnc.atleast(lits=list(ind), bound=k, top_id=nv + k,
                          encoding=EncType.seqcounter)
    sel0 = max([nv + k] + [abs(x) for cl in enc.clauses for x in cl]) + 1
    extra = list(enc.clauses)
    for i, (D, mem) in enumerate(rings):
        for c in range(k):
            extra.append([-(sel0 + i), -ind[c]] + [1 + v * k + c for v in mem])
    sv = Solver(name="cd15", bootstrap_with=cls + extra)
    capped = []
    for i, (D, mem) in enumerate(rings):
        if not sv.solve(assumptions=[sel0 + i]):
            capped.append((D, len(mem)))
    sv.delete()
    dens = len(E) / n
    print(f"  {label}: {n} pts, {len(E)} edges ({dens:.2f}/v), colours, "
          f"{len(rings)} rings, {len(capped)} CAPPED {capped[:5]}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    return capped


study(list(Sa), "Sa alone")
U, seen = list(Sa), set(Sa)
for D, name, rot in rots:
    fresh = [q for q in (rot(p) for p in Sa) if q not in seen]
    seen.update(fresh)
    U.extend(fresh)
    res = study(U, f"+ D={D}{name}")
    if res == "WIN" or res is None:
        break
print("DONE", flush=True)
