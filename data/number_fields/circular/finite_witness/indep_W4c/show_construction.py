#!/usr/bin/env python3
"""Print, for transparency, every constructed vector of claim 1 / claim 2
(exact, as integer numerators over D) next to the generator it matches and
the sign.  Usage: python3 show_construction.py WITNESS.json.gz"""
import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from biquad import Field  # noqa: E402
from check_geom import expected_q3_11, expected_q2_3  # noqa: E402

with gzip.open(sys.argv[1], "rt") as fh:
    d = json.load(fh)
D = d["denominator"]
if d["field"] == "Q(sqrt3, sqrt11)":
    ev, info = expected_q3_11(Field(3, 11), D)
    lab = "(j,k,l)"
else:
    ev, info = expected_q2_3(Field(2, 3), D)
    lab = "(j,l)"
gidx = {}
for i, g in enumerate(d["generators"]):
    gidx[tuple(g)] = (i, "+")
    gidx[tuple(-x for x in g)] = (i, "-")
used = set()
for key, v in sorted(ev.items()):
    m = gidx.get(v)
    if m:
        used.add(m[0])
    print("%s=%-12s vector*%d = %-40s matches generator #%s with sign %s" % (lab, key, D, list(v), m[0] if m else None, m[1] if m else None))
print("generators matched: %d of %d; constructed vectors: %d; unmatched constructed: %d"
      % (len(used), len(d["generators"]), len(ev), sum(1 for v in ev.values() if v not in gidx)))
print("info:", info)
