"""Bite the witness -- the operation that SHARPENS, which was never applied.

The spindle fails on the witness and the reason is measured: its disjunction
has strength exactly one, so each rotated copy needs a single coincidence and
may choose which.  But the spindle is stage TWO of de Grey's template, and
stage two consumes a NAMED pair.  What turns "one of three" into "this one" is
stage one and a half: the BITE.

They are different operations.  The spindle takes one copy per pair, rotated
about that pair's own end.  The bite takes ONE copy, rotated about a common
centre by the angle that makes a whole ring touch its own image -- Sa's
"one of three" becomes Y's named pair through rho_4 about the origin, and
nothing about that is per-pair.

So: for each centre worth trying and each ring the witness has about it, take
the biting rotation, union, and ask the exact forced-pair question of the
result.  A bite that produces a named forced pair at a closable distance is
de Grey's Y one level up, and the ordinary two-copy spindle finishes it.

Centres tried: every vertex of the witness, and the origin.  Rings: every
rational squared radius about that centre whose bite the field admits.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.field import Field
from hn.geometry import Point, rotation_joining, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance, forced_same
from pysat.solvers import Solver

k = 5
t0 = time.time()
doc = json.load(open("/home/user/darwin-50/research/hadwiger-nelson/"
                     "data/witness_five.json"))
K = Field(tuple(doc["field"]["generators"]))
W = [Point(K.element([Fr(c) for c in x]), K.element([Fr(c) for c in y]))
     for x, y in doc["points"]]
print(f"witness: {len(W)} points  [{time.time()-t0:.0f}s]", flush=True)


def rings(pts, C):
    out = defaultdict(int)
    for p in pts:
        d = p.dist2(C)
        if all(x == 0 for x in d.c[1:]) and d.c[0] != 0:
            out[d.c[0]] += 1
    return out


def study(U):
    basis = IntBasis.covering(U)
    rows = basis.rows(U)
    dim, D2 = basis.dim, basis.D * basis.D
    n = len(U)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    Eset = set(E)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return n, len(E), None
    pairs = []
    for i in range(n - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) in Eset:
                continue
            dd = Fr(int(sq[off, 0]), D2)
            if dd and closable_distance(dd):
                pairs.append((i, j, dd))
    forced = [(i, j, d) for i, j, d in pairs if forced_same(sv, i, j, k)]
    sv.delete()
    return n, len(E), (len(pairs), forced)


centres = [(f"vertex {i}", W[i]) for i in range(len(W))]
centres.append(("origin", Point(K.zero(), K.zero())))
tried = 0
for label, C in centres:
    R = rings(W, C)
    usable = [D for D in R if closable_distance(D) and D != 1 and R[D] >= 2]
    for D in sorted(usable, key=lambda d: -R[d])[:3]:
        rot = rotation_joining(D, K).about(C)
        inv = Rotation(rotation_joining(D, K).cos,
                       -rotation_joining(D, K).sin).about(C)
        seen, U = set(W), list(W)
        for r in (rot, inv):
            for p in W:
                q = r(p)
                if q not in seen:
                    seen.add(q)
                    U.append(q)
        n, m, res = study(U)
        tried += 1
        if res is None:
            print(f"*** {label}, ring D={D} ({R[D]} points): union {n} pts "
                  f"NOT 5-COLOURABLE -- chi >= 6 ***", flush=True)
            sys.exit(0)
        npairs, forced = res
        if forced:
            print(f"*** {label}, ring D={D}: {n} pts, {len(forced)} FORCED "
                  f"pairs: {forced[:3]} ***  [{time.time()-t0:.0f}s]",
                  flush=True)
            sys.exit(0)
        if tried % 15 == 0:
            print(f"   {tried} bites tried, latest {label} D={D} "
                  f"({R[D]} on the ring): {n} pts, {npairs} closable pairs, "
                  f"0 forced  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} bites, none forced a pair  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
