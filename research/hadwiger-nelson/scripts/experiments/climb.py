"""Run the same recipe at five colours, and skip the middle step.

The template that rebuilds de Grey's graph from Sa is:

    Z = C u rho(C)                       rho = rotation_joining(D) about the
                                         core's centre
    W = Z u sigma(Z)                     sigma = rotation_joining(4D) about
                                         one end of the ring's antipodal pair

and W is not k-colourable EXACTLY WHEN Z forces that pair.  So the forcing
never has to be proved on its own: build W and ask it once.  If Z forces, W is
uncolourable and that is the answer; if it does not, W colours and that is a
clean negative.  One SAT call replaces a forced-pair search over a union of
27000 points.

The core here is G*, the dihedral closure of G about its own pivot: 13873
points, 5-chromatic, exactly invariant -- Sa's role one level up.  Its two
doubly-usable rings are D = 4, de Grey's own, and D = 17/2, which Sa has no
analogue of.  W runs to four copies of G*, so about 55000 points.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import doubly_usable_ring
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)
G = build_G(F, as_graph=False)
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
print(f"G*: {len(Gs)} points  [{time.time()-t0:.0f}s]", flush=True)

for D in (Fr(4), Fr(17, 2)):
    assert doubly_usable_ring(D)
    rho = rotation_joining(D, F).about(PIV)
    Z = list(Gs)
    zs = set(Gs)
    for p in Gs:
        q = rho(p)
        if q not in zs:
            zs.add(q)
            Z.append(q)
    # one end of the ring's antipodal pair, as the spindle pivot
    ends = [p for p in Z
            if (p.x - PIV.x) ** 2 + (p.y - PIV.y) ** 2 == F.rational(D)]
    if not ends:
        print(f"  D={D}: the ring is empty in Z, skipping", flush=True)
        continue
    A = ends[0]
    B = Point(PIV.x + PIV.x - A.x, PIV.y + PIV.y - A.y)
    if B not in zs:
        print(f"  D={D}: no antipode in Z, skipping", flush=True)
        continue
    sigma = rotation_joining(4 * D, F).about(B)
    d2 = (A.x - sigma(A).x) ** 2 + (A.y - sigma(A).y) ** 2
    W = list(Z)
    ws = set(Z)
    for p in Z:
        q = sigma(p)
        if q not in ws:
            ws.add(q)
            W.append(q)
    print(f"  D={D}: Z has {len(Z)}, W has {len(W)} points; A and sigma(A) "
          f"at squared distance {d2}  [{time.time()-t0:.0f}s]", flush=True)
    g = build_graph(W)
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    n = len(W)
    print(f"  D={D}: {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)
    k = 5
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    cf = sv.accum_stats().get("conflicts", 0)
    sv.delete()
    print(f"  D={D}: W is 5-colourable {ok}  ({cf} conflicts)  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("  *** W IS NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
        with open(SC + f"climb_six_{D.numerator}_{D.denominator}.pkl",
                  "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in W], E), fh)
        break
