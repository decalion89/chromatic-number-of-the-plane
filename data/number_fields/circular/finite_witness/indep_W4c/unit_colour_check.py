#!/usr/bin/env python3
"""Informational: is the stored colouring also proper on the unit-distance
non-edge pairs (i.e. on the full unit-distance graph of the point set)?
Exact integer arithmetic.  Usage: python3 unit_colour_check.py WITNESS.json.gz"""
import gzip, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biquad import Field, int_sq_norm_8
with gzip.open(sys.argv[1], "rt") as fh:
    d = json.load(fh)
F = Field(3, 11) if d["field"] == "Q(sqrt3, sqrt11)" else Field(2, 3)
D = d["denominator"]
P = [tuple(p) for p in d["points"]]
E = set((min(e), max(e)) for e in d["edges"])
col = d["colouring"]
n = len(P)
extra = []
for i in range(n):
    for j in range(i + 1, n):
        dd = tuple(b - a for a, b in zip(P[i], P[j]))
        if int_sq_norm_8(F, dd) == (D * D, 0, 0, 0) and (i, j) not in E:
            extra.append((i, j))
mono = [(i, j) for (i, j) in extra if col[i] == col[j]]
print("unit-distance non-edge pairs: %d; monochromatic under the stored colouring: %d" % (len(extra), len(mono)))
print("examples:", mono[:5])
