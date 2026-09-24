"""The torus gate: a PERIODIC colouring of the physical plane that is proper for the module's
unit directions only.  f : R^2 -> [5], constant on the cells of an h-grid of the torus
R^2 / (P Z)^2; for every unit vector v of the module and every cell C, every cell meeting C + v
must get another colour (conservative, so any SAT answer is a genuine colouring).  Then
c(x) = f(position of x) properly colours the whole unit-distance graph of the module -- and it
has nothing to do with any coset colouring.  A 5-colouring here would refute rigidity outright.
usage: torusgate.py <module.json> <P> <G (cells per side)> [K]
"""
import sys, json, math, time
from fractions import Fraction as Fr
from pysat.solvers import Solver
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
d = json.load(open(f"/home/user/darwin-50/research/hadwiger-nelson/data/{sys.argv[1]}"))
P_ = float(sys.argv[2]); G = int(sys.argv[3]); K = int(sys.argv[4]) if len(sys.argv) > 4 else 5
F = Field(tuple(d["field_generators"]))
pts = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(pts)
U = sorted({(round((pts[b] - pts[a]).fx, 12), round((pts[b] - pts[a]).fy, 12)) for a, b in g.edges()}
           | {(round((pts[a] - pts[b]).fx, 12), round((pts[a] - pts[b]).fy, 12)) for a, b in g.edges()})
h = P_ / G
HEX = __import__("os").environ.get("HEX", "1") == "1"
# lattice basis of the torus: square (P,0),(0,P) or hexagonal (P,0),(P/2, P*sqrt3/2); cells are parallelograms
import numpy as _np
Bm = _np.array([[P_, 0.0], [0.0, P_]]) if not HEX else _np.array([[P_, P_ / 2], [0.0, P_ * math.sqrt(3) / 2]])
Binv = _np.linalg.inv(Bm)
X = lambda i, j, c: 1 + ((i % G) * G + (j % G)) * K + c
cl = [[X(i, j, c) for c in range(K)] for i in range(G) for j in range(G)]
pairs = set()
for vx, vy in U:
    al, be = Binv @ _np.array([vx, vy])
    fx, fy = al * G, be * G
    dx = [math.floor(fx)] + ([math.floor(fx) + 1] if abs(fx - round(fx)) > 1e-9 else [])
    dy = [math.floor(fy)] + ([math.floor(fy) + 1] if abs(fy - round(fy)) > 1e-9 else [])
    for a in dx:
        for b in dy:
            pairs.add((a % G, b % G))
if (0, 0) in pairs:
    print(f"{sys.argv[1]}: P={P_}, G={G}: a unit vector maps a cell onto itself -- impossible at this period"); sys.exit()
for i in range(G):
    for j in range(G):
        for a, b in pairs:
            i2, j2 = (i + a) % G, (j + b) % G
            if (i, j) < (i2, j2):
                for c in range(K): cl.append([-X(i, j, c), -X(i2, j2, c)])
cl.append([X(0, 0, 0)])
t0 = time.time()
s = Solver(name="cadical195", bootstrap_with=cl)
r = s.solve()
print(f"{sys.argv[1]}: {len(U)} unit vectors; torus P={P_} grid {G}x{G} (h={h:.3f}), {len(pairs)} offsets; "
      f"{K}-colourable: {r}   [{time.time()-t0:.1f}s]", flush=True)
if r:
    m = s.get_model()
    f = [[next(c for c in range(K) if m[X(i, j, c) - 1] > 0) for j in range(G)] for i in range(G)]
    json.dump({"P": P_, "G": G, "K": K, "f": f}, open(f"torus_{sys.argv[1]}_{P_}_{G}.json", "w"))
