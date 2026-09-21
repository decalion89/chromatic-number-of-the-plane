"""Keep the disjunction AND the rings, instead of trading one for the other.

Minimising the witness to 63 points kept everything load-bearing for "some
pair is monochromatic" and threw away the ring structure, because ring
structure is not load-bearing for that.  Then the bite had nothing to hook:
rings of two to four points against Sa's six, and 177 bites forced nothing.

The two requirements do not have to be traded.  G still has its twelve-point
rings.  So build the union of the obstruction's support -- the vertices that
carry the disjunction -- with the points of a populated ring about some
centre, and bite THAT.  The disjunction survives because its support is all
there; the ring survives because it was added back on purpose.

Tried over every centre of G that has a doubly-usable ring with enough points
on it, taking the ring and its neighbourhood alongside the 63-point support,
and asking the exact forced-pair question of each bitten union.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, forced_same, doubly_usable_ring
from pysat.solvers import Solver

k = 5
t0 = time.time()
doc = json.load(open("/home/user/darwin-50/research/hadwiger-nelson/"
                     "data/witness_five.json"))
from hn.field import Field
KW = Field(tuple(doc["field"]["generators"]))
Wpts = {Point(KW.element([Fr(c) for c in x]), KW.element([Fr(c) for c in y]))
        for x, y in doc["points"]}
G = build_G(K, as_graph=False)
support = [i for i, p in enumerate(G) if p in Wpts]
print(f"witness support inside G: {len(support)} of {len(G)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
basis = IntBasis.covering(G)
rows = basis.rows(G)
dim, D2 = basis.dim, basis.D * basis.D
EG = sorted(set((min(a, b), max(a, b))
                for a, b in fast_edges_complete(basis, rows)))
adj = defaultdict(set)
for a, b in EG:
    adj[a].add(b)
    adj[b].add(a)


def rings_about(ci):
    d = rows - rows[ci]
    sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
    rat = np.ones(len(sq), dtype=bool)
    for m in range(1, dim):
        rat &= sq[:, m] == 0
    out = defaultdict(list)
    for off in np.nonzero(rat)[0]:
        v = int(sq[off, 0])
        if v:
            out[Fr(v, D2)].append(int(off))
    return out


def study(U):
    b = IntBasis.covering(U)
    r = b.rows(U)
    d2 = b.D * b.D
    dm = b.dim
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    Eset = set(E)
    n = len(U)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return n, len(E), None
    pairs = []
    for i in range(n - 1):
        d = r[i + 1:] - r[i]
        sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dm):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) in Eset:
                continue
            dd = Fr(int(sq[off, 0]), d2)
            if dd and closable_distance(dd):
                pairs.append((i, j, dd))
    forced = [(i, j, d) for i, j, d in pairs if forced_same(sv, i, j, k)]
    sv.delete()
    return n, len(E), (len(pairs), forced)


best = []
for ci in range(len(G)):
    R = rings_about(ci)
    for D, members in R.items():
        if len(members) >= 6 and doubly_usable_ring(D):
            best.append((len(members), ci, D, members))
best.sort(reverse=True)
print(f"{len(best)} (centre, doubly-usable ring) with 6+ points; biggest "
      f"{[(b[0], str(b[2])) for b in best[:5]]}  [{time.time()-t0:.0f}s]",
      flush=True)

tried = 0
for size, ci, D, members in best[:40]:
    keep = set(support) | set(members) | {ci}
    for v in list(members) + [ci]:
        keep |= adj[v]
    U0 = [G[i] for i in sorted(keep)]
    rot = rotation_joining(D, K)
    f, iv = rot.about(G[ci]), Rotation(rot.cos, -rot.sin).about(G[ci])
    seen, U = set(U0), list(U0)
    for r in (f, iv):
        for p in U0:
            q = r(p)
            if q not in seen:
                seen.add(q)
                U.append(q)
    n, m, res = study(U)
    tried += 1
    if res is None:
        print(f"*** centre {ci}, ring D={D} ({size} points): {n} pts NOT "
              f"5-COLOURABLE -- chi >= 6 ***", flush=True)
        break
    npairs, forced = res
    if forced:
        print(f"*** centre {ci}, ring D={D}: {n} pts, {len(forced)} FORCED: "
              f"{forced[:3]} ***  [{time.time()-t0:.0f}s]", flush=True)
        break
    if tried % 8 == 0:
        print(f"   {tried} tried, latest centre {ci} D={D} ({size} on the "
              f"ring): base {len(U0)} -> union {n} pts, {npairs} closable "
              f"pairs, 0 forced  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} bites with the rings kept  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
