"""Record the first 5-chromatic unit-distance graph in this project that has a
symmetry group of its own.

Sa is glued at all six vertices of one C6 orbit at once -- so the union is
invariant rather than symmetrised afterwards -- and the forced pair that
appears is spindled over its whole orbit, six pivots at once, which is again
invariant because g rot_w g^-1 = rot_{g(w)}.  The result carries C6 by
construction, where every earlier graph here had isometry group of order 1.
"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from hn.coloring import is_k_colorable

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    return out
seen, U = set(Sa), list(Sa)
for w in orbit(Sa[25]):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
print(f"carrier n={len(U)}", flush=True)
spin = rotation_joining(Fr(64, 9), F)
seen2, V = set(U), list(U)
for w in orbit(U[1]):
    rot = spin.about(w)
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
g = build_graph(V); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"symmetric spindle: n={n} m={m} deg={2.0*m/n:.2f}", flush=True)
t0 = time.time()
ok4, _ = is_k_colorable(g, 4, timeout=1800)
print(f"  4-colourable: {ok4}   [{time.time()-t0:.0f}s]", flush=True)
assert ok4 is False
json.dump({"field_generators": list(F.gens), "n": n, "m": m,
           "recipe": [
               "Sa = de Grey's 39 points under <rot60, reflect>, 397 points",
               "glued at ALL SIX vertices of the C6 orbit of Sa[25] at once,",
               "  each by the 60-degree rotation about that vertex -> 1021 points,",
               "  mean degree 13.34, rigid at four colours, 153 forced pairs",
               "spindled at d^2 = 64/9 over the whole C6 orbit of the pivot,",
               "  six rotations at once, so the result is C6-invariant",
           ],
           "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                       [[c.numerator, c.denominator] for c in p.y.c]]
                      for p in g.vertices]},
          open("/home/user/darwin-50/research/hadwiger-nelson/data/five_symmetric.json", "w"))
print("  written data/five_symmetric.json", flush=True)
