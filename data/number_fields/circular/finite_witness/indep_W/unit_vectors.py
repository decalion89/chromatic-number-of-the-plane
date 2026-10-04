#!/usr/bin/env python3
"""Count all unit vectors ((a + b r)/D, (c + e r)/D), r = sqrt d, by brute force (referee code), and compare with the
differences that occur as edges of the witness.  Usage: unit_vectors.py W.json.gz"""
import gzip, json, sys, math
W = json.load(gzip.open(sys.argv[1], 'rt')); d, D = W['d'], W['denominator']
B = math.isqrt(D * D // d)
U = set()
for b in range(-B, B + 1):
    for e in range(-B, B + 1):
        r = D * D - d * (b * b + e * e)
        if r < 0:
            continue
        for a in range(-D, D + 1):
            c2 = r - a * a
            if c2 < 0:
                continue
            c = math.isqrt(c2)
            if c * c == c2:
                for cc in {c, -c}:
                    if a * b + cc * e == 0:
                        U.add((a, b, cc, e))
P = W['points']
used = set()
for i, j in W['edges']:
    v = tuple(P[j][k] - P[i][k] for k in range(4)); used.add(v); used.add(tuple(-t for t in v))
rational = sum(1 for u in U if u[1] == 0 and u[3] == 0)
print(f"d={d} D={D}: {len(U)} unit vectors with denominator D ({rational} rational); {len(used)} occur as edge "
      f"differences; all edge differences among them: {used <= U}")
