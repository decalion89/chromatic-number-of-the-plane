#!/usr/bin/env python3
"""Informational: how many edges of H use each generator direction (+-g),
and the angle of each generator (float, informational only)."""
import gzip
import json
import math
import sys

with gzip.open(sys.argv[1], "rt") as fh:
    d = json.load(fh)
A, B = (3, 11) if "11" in d["field"] else (2, 3)
sa, sb = math.sqrt(A), math.sqrt(B)
D = d["denominator"]
P = [tuple(p) for p in d["points"]]
cnt = {}
for g in d["generators"]:
    cnt[tuple(g)] = 0
for (i, j) in d["edges"]:
    dd = tuple(b - a for a, b in zip(P[i], P[j]))
    nd = tuple(-x for x in dd)
    if dd in cnt:
        cnt[dd] += 1
    elif nd in cnt:
        cnt[nd] += 1
    else:
        print("edge with non-generator difference", i, j)


def fl(t):
    return (t[0] + t[1] * sa + t[2] * sb + t[3] * sa * sb) / D


unused = 0
for g, c in cnt.items():
    ang = math.degrees(math.atan2(fl(g[4:]), fl(g[:4])))
    print("generator %-40s angle %9.4f deg  edges %d" % (list(g), ang, c))
    if c == 0:
        unused += 1
print("generators used by no edge:", unused, "of", len(cnt))
