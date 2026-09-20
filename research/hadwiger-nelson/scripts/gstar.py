"""Is the symmetric closure of G still 5-colourable?

G is 5-chromatic and lopsided.  G* is its dihedral closure about (-2,0), the
pivot its two copies of Y were turned about -- 13873 points where twelve
disjoint copies would be 18972, so five thousand points are shared and the
copies genuinely interlock.  If G* were not 5-colourable the search would be
over.  If it is, then G* is a 5-chromatic graph that is exactly invariant
under the 12-element dihedral group, which is what Sa is one level down, and
rotating it about its own centre is the step that built Y out of Sa.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
r60 = _rot60(F)
PIV = Point(F.rational(-2), F.zero())
rot = r60.about(PIV)

G = build_G(F, as_graph=False)
out, seen = [], set()
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
                out.append(q)
print(f"G*: {len(out)} points  [{time.time()-t0:.0f}s]", flush=True)
g = build_graph(out)
E = [(min(a, b), max(a, b)) for a, b in g.edges()]
n = len(out)
print(f"G*: {n} points, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)
with open(SC + "gstar.pkl", "wb") as fh:
    pickle.dump(([(str(p.x), str(p.y)) for p in out], E), fh)

for k in (5, 4):
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    print(f"  G* is {k}-colourable: {ok}  "
          f"[{time.time()-t0:.0f}s, {sv.accum_stats().get('conflicts', 0)} "
          f"conflicts]", flush=True)
    sv.delete()
    if k == 5 and not ok:
        print("  *** G* IS NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
        break
