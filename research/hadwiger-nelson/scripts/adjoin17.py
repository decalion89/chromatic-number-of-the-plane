"""Go past where de Grey's chain stopped, by adjoining the root it wanted.

His doubling chain is D = 1, 3, 4, 16, and it halts at 16 because the next
step needs sqrt(16*16 - 1) = sqrt(255) = sqrt(3*5*17) and his field has no
sqrt(17).  The field was forced on him by the chain, not chosen; nothing
obliges a search to stay inside it.

The ring D = 16 of G* costs its two radicals separately.  rho needs
sqrt(4*16-1) = sqrt(63) = 3 sqrt(7), which is already there -- rho is 31/32 +
i (3/32) sqrt(7), the same turn that served as the spindle at D = 4.  Only the
spindle of the antipodal pair, which sits at squared distance 64, needs
sqrt(255).  So the extension is needed for one step at the very end, and the
copy it produces is still an exact unit-distance graph: its coordinates simply
live in Q(sqrt3, sqrt5, sqrt7, sqrt11, sqrt17), a 32-dimensional field instead
of a 16-dimensional one.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import (DEGREY_FIELD as K, Point, _rot60, rotation_joining)
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
K17 = Field((3, 5, 7, 11, 17))
RAD17 = (3, 5, 7, 11, 17)
t0 = time.time()

# Which rings pay both costs once sqrt(17) is available?
ALL = [Fr(1), Fr(4), Fr(31, 3), Fr(17, 2), Fr(7), Fr(20, 3), Fr(13, 3),
       Fr(7, 3), Fr(3, 2), Fr(16)]
for D in ALL:
    a = closable_distance(D, RAD17)
    b = closable_distance(4 * D, RAD17)
    if a and b and D != 1:
        print(f"  over K(sqrt17), D={D} pays both", flush=True)

PIV17 = Point(K17.rational(-2), K17.zero())
rot = _rot60(K).about(Point(K.rational(-2), K.zero()))
G = build_G(K, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
Gs = [Point(embed(p.x, K17), embed(p.y, K17)) for p in Gs]
print(f"G* embedded in K(sqrt17): {len(Gs)} points  "
      f"[{time.time()-t0:.0f}s]", flush=True)

D = Fr(16)
rho = rotation_joining(D, K17).about(PIV17)
Z, zs = list(Gs), set(Gs)
for p in Gs:
    q = rho(p)
    if q not in zs:
        zs.add(q)
        Z.append(q)
ring = [p for p in Z
        if (p.x - PIV17.x) ** 2 + (p.y - PIV17.y) ** 2 == K17.rational(D)]
anti = []
for p in ring:
    q = Point(PIV17.x + PIV17.x - p.x, PIV17.y + PIV17.y - p.y)
    if q in zs:
        anti.append((p, q))
print(f"Z = G* u rho_16(G*): {len(Z)} points, {len(ring)} on the ring, "
      f"{len(anti)} antipodal  [{time.time()-t0:.0f}s]", flush=True)
if not anti:
    print("  no antipodal pair on the ring -- nothing to spindle", flush=True)
    sys.exit(0)

A, B = anti[0]
sigma = rotation_joining(4 * D, K17).about(B)
d2 = (A.x - sigma(A).x) ** 2 + (A.y - sigma(A).y) ** 2
W, ws = list(Z), set(Z)
for p in Z:
    q = sigma(p)
    if q not in ws:
        ws.add(q)
        W.append(q)
print(f"W = Z u sigma_64(Z): {len(W)} points; A and sigma(A) at squared "
      f"distance {d2}  [{time.time()-t0:.0f}s]", flush=True)
g = build_graph(W)
E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
n = len(W)
print(f"W: {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)
k = 5
cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
for a, b in E:
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
ok = sv.solve()
cf = sv.accum_stats().get("conflicts", 0)
sv.delete()
print(f"W is 5-colourable: {ok}  ({cf} conflicts)  [{time.time()-t0:.0f}s]",
      flush=True)
if not ok:
    print("  *** NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
    with open(SC + "adjoin17_six.pkl", "wb") as fh:
        pickle.dump(([(str(p.x), str(p.y)) for p in W], E), fh)
