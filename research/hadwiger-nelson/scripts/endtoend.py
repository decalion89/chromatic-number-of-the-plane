"""Run the whole recipe one level down, from the criterion, and check it lands.

Nothing here is taken from de Grey's answer.  The ring is chosen by the
criterion -- of Sa's five rational closable rings, D = 4 is the only one that
pays both costs -- and everything after follows from it mechanically:

    core     Sa, the dihedral orbit of S: 397 points, invariant
    rho      rotation_joining(4) about the centre: the turn that makes the
             ring of squared radius 4 produce cross edges
    Z        Sa u rho(Sa), which should force the ring's antipodal pair
    pivot    one end of that pair, (-2,0)
    sigma    rotation_joining(16) about the pivot: the pair is at squared
             distance 16, so this is the turn that brings the FAR end of the
             pair to distance 1 from itself
    W        Z u sigma(Z)

If the pair is forced in Z then in W the far end A agrees with the pivot B,
and sigma(A) agrees with B too, because sigma fixes B and sigma(Z) is a copy
of Z.  So A and sigma(A) agree -- and they are at distance 1.  W cannot be
4-coloured.

That is a proof, not a search, and the only input it needs about Sa is that
its ring pair is forced.  If W comes back non-4-colourable, the pipeline
reproduces chi >= 5 end to end from the criterion alone, and the same template
is what gets applied at five colours.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import doubly_usable_ring
from pysat.solvers import Solver

t0 = time.time()
D = Fr(4)
assert doubly_usable_ring(D), "the criterion must pick this ring"

Sa = build_Sa(F)
rho = rotation_joining(D, F)          # about the origin: Sa's own centre
Z, seen = list(Sa), set(Sa)
for p in Sa:
    q = rho(p)
    if q not in seen:
        seen.add(q)
        Z.append(q)
print(f"Z = Sa u rho(Sa): {len(Z)} points  [{time.time()-t0:.0f}s]",
      flush=True)

A = Point(F.rational(2), F.zero())
B = Point(F.rational(-2), F.zero())
assert A in seen and B in seen, "the antipodal pair must be present"

sigma = rotation_joining(4 * D, F).about(B)   # 4D = 16: the pair's distance
assert sigma(B) == B, "the pivot must be fixed"
W, seen2 = list(Z), set(Z)
for p in Z:
    q = sigma(p)
    if q not in seen2:
        seen2.add(q)
        W.append(q)
g = build_graph(W)
E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
n = len(W)
iA, iAp = W.index(A), W.index(sigma(A))
d2 = ((W[iA].x - W[iAp].x) ** 2 + (W[iA].y - W[iAp].y) ** 2)
print(f"W = Z u sigma(Z): {n} points, {len(E)} edges; "
      f"A and sigma(A) are at squared distance {d2}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

for k in (4, 3):
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    cf = sv.accum_stats().get("conflicts", 0)
    sv.delete()
    print(f"  W is {k}-colourable: {ok}  ({cf} conflicts) "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if k == 4 and not ok:
        print("  *** chi(W) >= 5 -- the recipe lands, from the criterion "
              "alone ***", flush=True)
        break
