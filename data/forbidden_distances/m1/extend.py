"""exact extension of a configuration by candidate points (given as complex floats); checks every new point in
Q(sqrt3, sqrt11) and at distance 1 from at least two points. usage: python3 extend.py BASE.json OUT.json re,im ..."""
import sys
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exactconf as ec
pts, _ = ec.load(sys.argv[1])
for s in sys.argv[3:]:
    re, im = (float(v) for v in s.split(','))
    p = ec.refine(pts, complex(re, im), tol=1e-3)
    pts.append(p)
    d2, edges = ec.graph(pts)
    print('added', s, '->', ec.kstr(p[0]), ec.kstr(p[1]), 'unit degree', sum(1 for e in edges if len(pts) - 1 in e))
ec.save(sys.argv[2], pts)
d2, edges = ec.graph(pts); atoms = ec.independent_sets(len(pts), edges)
print(len(pts), 'points', len(edges), 'unit pairs', len(atoms), 'independent sets')
