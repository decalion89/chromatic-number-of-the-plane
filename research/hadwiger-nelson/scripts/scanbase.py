"""Which base vertices are worth gluing at?

The symmetric glue is parametrised by which orbit of centres it uses, and the
orbits are not equal: gluing at the origin is a symmetry of Sa and buys
nothing, while the orbit of Sa[25] takes the mean degree from 9.94 to 13.34.
This ranks every base vertex by the density its own orbit produces, with no
solver in the loop, so the expensive scans can be pointed at the good ones.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out

seenorb = set()
rows = []
t0 = time.time()
for b in range(len(Sa)):
    O = orbit(Sa[b])
    key = frozenset(O)
    if key in seenorb: continue
    seenorb.add(key)
    seen, U = set(Sa), list(Sa)
    for w in O:
        rot = g60.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    if len(U) == len(Sa): continue
    g = build_graph(U); m = sum(len(a) for a in g.adj) // 2
    rows.append((2.0 * m / g.n, g.n, m, b, len(O)))
rows.sort(reverse=True)
print(f"{len(rows)} distinct orbits of centres   [{time.time()-t0:.0f}s]", flush=True)
print(f"{'deg':>6} {'n':>6} {'m':>7} {'base':>6} {'centres':>8}")
for r in rows[:15]:
    print(f"{r[0]:>6.2f} {r[1]:>6} {r[2]:>7} v{r[3]:<5} {r[4]:>8}")
print("  ...")
for r in rows[-3:]:
    print(f"{r[0]:>6.2f} {r[1]:>6} {r[2]:>7} v{r[3]:<5} {r[4]:>8}")
