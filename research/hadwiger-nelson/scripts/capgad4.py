"""Peel Sa by ORBIT, keeping the symmetry that makes the cap work.

The cap -- Sa's D = 4 ring takes at most two colours in every 4-colouring --
is an UNSAT proof and costs about four minutes however it is asked, so peeling
one vertex at a time would take a day and the unsatisfiable core comes back
naming all 397 vertices, which discriminates nothing.

But Sa is the twelve-element dihedral closure of a 39-point seed, and the cap
is a statement about a ring that the same group preserves.  So the natural
unit is the ORBIT, not the vertex: thirty-odd trials per sweep instead of four
hundred, and every intermediate stays as symmetric as Sa is, which is the
property the construction was built around in the first place.

What survives is the symmetric core of de Grey's gadget -- the thing worth
trying to rebuild one level up.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Point, _rot60
from hn.fast import IntBasis, fast_edges_complete
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 2
t0 = time.time()
P = build_Sa(K)
n = len(P)
idx = {p: i for i, p in enumerate(P)}
rot60 = _rot60(K)
ZERO = Point(K.zero(), K.zero())

# dihedral orbits: rotate six ways, reflect in y
orb = {}
orbits = []
for i, p in enumerate(P):
    if i in orb:
        continue
    members = set()
    for base in (p, Point(p.x, -p.y)):
        q = base
        for _ in range(6):
            members.add(idx[q])
            q = rot60(q)
    o = len(orbits)
    orbits.append(sorted(members))
    for m in members:
        orb[m] = o
sizes = defaultdict(int)
for o in orbits:
    sizes[len(o)] += 1
print(f"Sa: {n} points in {len(orbits)} dihedral orbits, sizes "
      f"{dict(sorted(sizes.items()))}  [{time.time()-t0:.0f}s]", flush=True)

b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
dm, d2 = b.dim, b.D * b.D
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
zi = idx[ZERO]
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
print(f"ring {ring} lives in orbit(s) {sorted({orb[i] for i in ring})}",
      flush=True)

nv = n * k
sel = [nv + 1 + v for v in range(n)]
ind = [nv + n + 1 + c for c in range(k)]
top = nv + n + k
cls = []
for v in range(n):
    cls.append([-sel[v]] + [1 + v * k + c for c in range(k)])
    for c in range(k):
        cls.append([sel[v], -(1 + v * k + c)])
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
for c in range(k):
    cls.append([-ind[c]] + [1 + v * k + c for v in ring])
enc = CardEnc.atleast(lits=list(ind), bound=CAP + 1, top_id=top,
                      encoding=EncType.seqcounter)
sv = Solver(name="cd15", bootstrap_with=cls + list(enc.clauses))
print(f"solver built  [{time.time()-t0:.0f}s]", flush=True)

pinned = {orb[zi]} | {orb[i] for i in ring}


def holds(keep_orbits):
    pres = {v for o in keep_orbits for v in orbits[o]}
    a = [sel[v] if v in pres else -sel[v] for v in range(n)]
    return not sv.solve(assumptions=a)


present = set(range(len(orbits)))
assert holds(present), "the cap does not hold on Sa itself"
print(f"cap confirmed  [{time.time()-t0:.0f}s]", flush=True)
rounds = 0
while True:
    rounds += 1
    dropped = 0
    for o in sorted(present):
        if o in pinned or o not in present:
            continue
        if holds(present - {o}):
            present.discard(o)
            dropped += 1
            pts = sum(len(orbits[x]) for x in present)
            print(f"    dropped orbit {o} ({len(orbits[o])} pts) -> {pts} "
                  f"points  [{time.time()-t0:.0f}s]", flush=True)
    pts = sum(len(orbits[x]) for x in present)
    print(f"  sweep {rounds}: dropped {dropped} orbits, now {len(present)} "
          f"orbits / {pts} points  [{time.time()-t0:.0f}s]", flush=True)
    if not dropped:
        break
keep = sorted(v for o in present for v in orbits[o])
json.dump({"k": k, "cap": CAP, "orbits": len(present),
           "points": [[[str(x) for x in P[i].x.c],
                       [str(y) for y in P[i].y.c]] for i in keep]},
          open("capgadget.json", "w"))
print(f"\nsymmetric core: {len(keep)} points in {len(present)} orbits"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
