"""Rebuild Exoo-Ismailescu's 6-chromatic {1,2}-graph from the paper and re-verify it.

arXiv 1909.13177 (Geombinatorics 2020): every 5-colouring of the plane has two
points of the same colour at distance 1 or 2.  Vertices are
[a,b,c,d] = (a sqrt3/12 + b sqrt11/12,  c/12 + d sqrt33/12), i.e. our field.
Claims to check: G 205 vertices / 966 unit / 423 two-edges; H 214 / 1004 / 446
with A,B (distance 5) the same colour in every 5-colouring; K = H u rot_A(H)
(cos 49/50) 426 / 2009 / 892 and not 5-colourable.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json, itertools, time
from fractions import Fraction as Fr
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point
from pysat.solvers import Solver
t0 = time.time()
F = Field((3, 11)); r3, r11 = F.sqrt(3), F.sqrt(11); r33 = r3 * r11
q = lambda v: F.rational(Fr(v, 12))
def P(a, b, c, d): return Point(q(a) * r3 + q(b) * r11, q(c) + q(d) * r33)
S = [(0,0,0,0),(0,0,0,-4),(0,0,-6,-2),(0,0,-6,2),(-6,0,0,-2),(-4,0,0,0),(-4,0,-6,-2),
     (-4,0,-6,2),(-2,0,0,-2),(-2,0,-6,-4),(-2,0,-6,4),(0,-6,-6,0),(-5,-3,3,3),(-5,3,-3,3),
     (-2,-6,0,0),(-2,-6,0,-4),(-2,6,0,0),(-2,6,0,-4),(-6,-6,0,0),(-6,6,0,0),(-4,0,0,-4),
     (0,0,-12,0),(-8,0,0,0)]
assert len(S) == 23
T = set()
for a, b, c, d in S:
    T.add((a, b, c, d)); T.add((a, b, -c, -d)); T.add((-a, -b, c, d))
print("T:", len(T), "(paper: 57)")
def rot(p):  # the paper's formula for the pi/3 rotation
    a, b, c, d = p
    return ((a - c) // 2 if (a - c) % 2 == 0 else Fr(a - c, 2), (b - 3 * d) / 2,
            (3 * a + c) / 2, (b + d) / 2)
half = F.rational(Fr(1, 2)); c60, s60 = half, r3 * half
def R(p): return Point(p.x * c60 - p.y * s60, p.x * s60 + p.y * c60)
V = {}
for p in T:
    pt = P(*p)
    for k in range(6):
        V[pt] = True; pt = R(pt)
V = list(V)
def graph(pts):
    one, four = F.rational(1), F.rational(4)
    E1, E2 = [], []
    for i, j in itertools.combinations(range(len(pts)), 2):
        dx = pts[i].fx - pts[j].fx; dy = pts[i].fy - pts[j].fy; dd = dx*dx + dy*dy
        if abs(dd - 1) < 1e-9 and pts[i].dist2(pts[j]) == one: E1.append((i, j))
        elif abs(dd - 4) < 1e-9 and pts[i].dist2(pts[j]) == four: E2.append((i, j))
    return E1, E2
E1, E2 = graph(V)
print(f"G: {len(V)} vertices, {len(E1)} unit, {len(E2)} two-edges   (paper 205/966/423)")
extra = [(-2,0,0,-6),(8,0,0,4),(-4,-6,-6,-4),(-4,6,6,-4),(-3,-3,-3,-5),(-4,0,-12,4),
         (-4,0,12,4),(7,-3,3,3),(7,3,-3,3)]
A, B = P(*extra[0]), P(*extra[1])
print("AB^2 =", A.dist2(B))
VH = list(dict.fromkeys(V + [P(*e) for e in extra]))
E1h, E2h = graph(VH)
print(f"H: {len(VH)} vertices, {len(E1h)} unit, {len(E2h)} two-edges   (paper 214/1004/446)")
def cnf_of(n, edges, K=5):
    # symmetry breaking: a triangle of the graph gets colours 0, 1, 2 (sound: any
    # proper colouring can be permuted so), found among the edges given
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in edges:
        for c in range(K): cl.append([-X(a, c), -X(b, c)])
    nb = {}
    for a, b in edges: nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
    tri = next(((a, b, c) for a, b in edges for c in nb[a] & nb[b]), None)
    if tri:
        for k, v in enumerate(tri): cl.append([X(v, k)])
    return cl, X
cl, X = cnf_of(len(VH), E1h + E2h)
ia, ib = VH.index(A), VH.index(B)
s = Solver(name="cd19", bootstrap_with=cl); print("H 5-colourable:", s.solve())
for c in range(5): s.add_clause([-X(ia, c), -X(ib, c)])
print("H with A != B 5-colourable:", s.solve(), "(paper: False, A=B forced)")
ca, sa = F.rational(Fr(49, 50)), F.rational(Fr(3, 50)) * r11
def RA(p):
    dx, dy = p.x - A.x, p.y - A.y
    return Point(A.x + dx * ca - dy * sa, A.y + dx * sa + dy * ca)
Bp = RA(B); print("BB'^2 =", B.dist2(Bp))
VK = list(dict.fromkeys(VH + [RA(p) for p in VH]))
E1k, E2k = graph(VK)
print(f"K: {len(VK)} vertices, {len(E1k)} unit, {len(E2k)} two-edges   (paper 426/2009/892)")
for name in ("cd19", "g4", "m22"):
    cl, X = cnf_of(len(VK), E1k + E2k)
    s = Solver(name=name, bootstrap_with=cl)
    print(f"K 5-colourable ({name}):", s.solve(), f"  [{time.time()-t0:.0f}s]", flush=True)
cl, X = cnf_of(len(VK), E1k)
s = Solver(name="cd19", bootstrap_with=cl)
print("K with unit edges only 5-colourable:", s.solve())
OUT = os.environ.get("HN_OUT", "/tmp/hn")
os.makedirs(OUT, exist_ok=True)
json.dump({"field_generators": [3, 11],
           "points": [[[[t.numerator, t.denominator] for t in p.x.c],
                       [[t.numerator, t.denominator] for t in p.y.c]] for p in VK],
           "unit_edges": E1k, "two_edges": E2k},
          open(os.path.join(OUT, "exoo_ismailescu_K426.json"), "w"))
print("written", os.path.join(OUT, "exoo_ismailescu_K426.json"))
