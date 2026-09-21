"""Write the witness out, and verify it from the file alone.

The obstruction is now small enough to hand over: sixty-odd vertices of G,
their unit edges, and the forbidden pairs, such that no proper 5-colouring
avoids all of the latter.  Stated that way it is a claim about a search
nobody else ran.  Written to a file with exact coordinates and re-checked
from that file, it is a claim anyone can test.

So: regenerate the witness, dump it as JSON with each coordinate as its exact
vector over the field basis, and then verify -- reading only the file --
that every listed edge really is at distance one, every forbidden pair really
is at the distance it claims, and the graph with those pairs added is not
5-colourable.  The verification uses the file's numbers, not the search's
memory.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from hn.degrey import build_G
from hn.field import Field
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 5
t0 = time.time()
OUT = "/home/user/darwin-50/research/hadwiger-nelson/data/witness_five.json"
KEEP = [Fr(15, 16), Fr(16), Fr(17, 2), Fr(9), Fr(7)]
P = build_G(K, as_graph=False)
basis = IntBasis.covering(P)
rows = basis.rows(P)
dim, D2 = basis.dim, basis.D * basis.D
n = len(P)
E = sorted(set((min(a, b), max(a, b))
               for a, b in fast_edges_complete(basis, rows)))
Eset = set(E)
byd = defaultdict(list)
for i in range(n - 1):
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
                byd[Fr(v, D2)].append((i, j))
FORB = [(D, p) for D in KEEP for p in byd[D]]
verts = sorted({v for _, (a, b) in FORB for v in (a, b)})


def colours(vs, prs):
    vs = set(vs)
    ren = {v: i for i, v in enumerate(sorted(vs))}
    cls = [[1 + i * k + c for c in range(k)] for i in range(len(ren))]
    for a, b in E:
        if a in vs and b in vs:
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    for _, (a, b) in prs:
        if a in vs and b in vs:
            for c in range(k):
                cls.append([-(1 + ren[a] * k + c), -(1 + ren[b] * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


cur = set(verts)
while True:
    dropped = 0
    for v in sorted(cur):
        if v not in cur:
            continue
        if not colours(cur - {v}, FORB):
            cur -= {v}
            dropped += 1
    if dropped == 0:
        break
FORB = [x for x in FORB if x[1][0] in cur and x[1][1] in cur]
while True:
    dropped = 0
    for target in sorted(FORB, key=lambda x: -float(x[0])):
        if target not in FORB:
            continue
        trial = [x for x in FORB if x is not target]
        if colours(cur, trial):
            continue
        FORB = trial
        dropped += 1
    if dropped == 0:
        break
verts = sorted(cur)
print(f"witness: {len(verts)} vertices, {len(FORB)} forbidden pairs, "
      f"colours {colours(verts, FORB)}  [{time.time()-t0:.0f}s]", flush=True)

ren = {v: i for i, v in enumerate(verts)}
sub_edges = [[ren[a], ren[b]] for a, b in E if a in cur and b in cur]
doc = {
    "field": {"generators": [3, 5, 7, 11],
              "basis": "squarefree products of the generators, in the order "
                       "hn.field.Field((3,5,7,11)) lists them"},
    "claim": "no proper 5-colouring of the unit-distance graph on these "
             "points leaves every listed pair bichromatic",
    "points": [[[str(c) for c in P[v].x.c], [str(c) for c in P[v].y.c]]
               for v in verts],
    "unit_edges": sub_edges,
    "forbidden_pairs": [[ren[a], ren[b], str(D)] for D, (a, b) in FORB],
    "classes": {str(D): sum(1 for DD, _ in FORB if DD == D) for D in KEEP},
}
import os
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(doc, fh, indent=1)
print(f"written to {OUT}  [{time.time()-t0:.0f}s]", flush=True)

# --- verify from the file alone -------------------------------------------
doc = json.load(open(OUT))
F2 = Field(tuple(doc["field"]["generators"]))
pts = [Point(F2.element([Fr(c) for c in x]), F2.element([Fr(c) for c in y]))
       for x, y in doc["points"]]
m = len(pts)
bad_e = sum(1 for a, b in doc["unit_edges"] if pts[a].dist2(pts[b]) != 1)
bad_p = sum(1 for a, b, D in doc["forbidden_pairs"]
            if pts[a].dist2(pts[b]) != Fr(D))
notclo = sum(1 for _, _, D in doc["forbidden_pairs"]
             if not closable_distance(Fr(D)))
cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
for a, b in doc["unit_edges"]:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
for a, b, _ in doc["forbidden_pairs"]:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="cd15", bootstrap_with=cls)
ok = sv.solve()
sv.delete()
print(f"\nfrom the file alone: {m} points, {len(doc['unit_edges'])} unit "
      f"edges, {len(doc['forbidden_pairs'])} forbidden pairs", flush=True)
print(f"   edges at the wrong distance: {bad_e}", flush=True)
print(f"   pairs at the wrong distance: {bad_p}", flush=True)
print(f"   pairs at a non-closable distance: {notclo}", flush=True)
print(f"   5-colourable: {ok}", flush=True)
print(f"   {'VERIFIED' if not (bad_e or bad_p or notclo or ok) else 'FAILS'}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
