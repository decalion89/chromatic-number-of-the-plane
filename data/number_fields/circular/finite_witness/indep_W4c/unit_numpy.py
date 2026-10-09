#!/usr/bin/env python3
"""Second, independent implementation of the exact unit-distance count,
vectorised with numpy int64.  Exactness: we bound every intermediate
quantity by an explicit a-priori bound and assert it is < 2^62, so no
int64 overflow can happen and all arithmetic is exact integer arithmetic.

Usage: python3 unit_numpy.py WITNESS.json.gz
"""
import gzip
import json
import sys

import numpy as np

with gzip.open(sys.argv[1], "rt") as fh:
    d = json.load(fh)
A, B = (3, 11) if d["field"] == "Q(sqrt3, sqrt11)" else (2, 3)
assert d["field"] in ("Q(sqrt3, sqrt11)", "Q(sqrt2, sqrt3)")
D = d["denominator"]
P = np.array(d["points"], dtype=np.int64)
n = len(P)
M = int(np.abs(P).max())
dm = 2 * M                          # bound on |difference coordinate|
bound = 2 * (1 + A + B + A * B) * dm * dm   # bound on every norm component
assert bound < 2 ** 62, bound
print("max |coord| = %d, a-priori bound on norm components = %d (< 2^62)" % (M, bound))

G = set()
for g in d["generators"]:
    G.add(tuple(g))
    G.add(tuple(-x for x in g))
E = set((min(e), max(e)) for e in d["edges"])

unit_pairs = []
for i in range(n - 1):
    Q = P[i + 1:] - P[i]
    a0, a1, a2, a3, b0, b1, b2, b3 = (Q[:, k] for k in range(8))
    c0 = a0 * a0 + A * a1 * a1 + B * a2 * a2 + A * B * a3 * a3 + b0 * b0 + A * b1 * b1 + B * b2 * b2 + A * B * b3 * b3
    c1 = 2 * (a0 * a1 + B * a2 * a3 + b0 * b1 + B * b2 * b3)
    c2 = 2 * (a0 * a2 + A * a1 * a3 + b0 * b2 + A * b1 * b3)
    c3 = 2 * (a0 * a3 + a1 * a2 + b0 * b3 + b1 * b2)
    hit = np.nonzero((c0 == D * D) & (c1 == 0) & (c2 == 0) & (c3 == 0))[0]
    for h in hit:
        unit_pairs.append((i, i + 1 + int(h)))
U = set(unit_pairs)
gen_unit = sum(1 for (i, j) in U if tuple(int(x) for x in (P[j] - P[i])) in G)
print("unit-distance pairs: %d; with generator difference: %d; listed edges among them: %d; unit non-edges: %d; edges not unit: %d"
      % (len(U), gen_unit, len(U & E), len(U - E), len(E - U)))
