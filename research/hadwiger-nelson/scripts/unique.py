"""Is the colouring rigid?  Unique k-colourability forces chi_c = k.

Sa and Y measure tight at 4 and the Moser spindle does not, and there is a
classical reason a graph can be tight: a uniquely k-colourable graph has
chi_c = k.  If Sa's 4-colouring is rigid -- if the partition into classes is
forced, not merely one of many -- then its tightness is explained rather
than observed, and the target one level up becomes concrete: a uniquely
5-colourable unit-distance graph would have chi_c = 5.

Rigidity is tested pair by pair with the forcing machinery already here.  In
a uniquely colourable graph every same-class pair is forced to agree and
every cross-class pair is forced to differ, so sampling both kinds measures
how close the colouring is to rigid.
"""
import sys, time, random, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
WHICH = sys.argv[1]
KC = int(sys.argv[2])
SAMP = int(sys.argv[3]) if len(sys.argv) > 3 else 120
t0 = time.time()
P = {"Sa": build_Sa, "Y": build_Y,
     "G": lambda k: build_G(k, as_graph=False)}[WHICH](K)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n = len(P)
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, c in E:
    for col in range(KC):
        cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve(), f"{WHICH} is not {KC}-colourable"
m = s.get_model()
col = {}
for v in range(n):
    for c in range(KC):
        if m[v * KC + c] > 0:
            col[v] = c
            break
sizes = [sum(1 for v in col if col[v] == c) for c in range(KC)]
print(f"{WHICH}: {n} pts, {len(E)} edges, {KC}-colouring classes {sizes}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
random.seed(11)
same = [(u, v) for u in range(n) for v in range(u + 1, n) if col[u] == col[v]]
diff = [(u, v) for u in range(n) for v in range(u + 1, n) if col[u] != col[v]]
random.shuffle(same)
random.shuffle(diff)


# forced_different(i, j): no colouring gives both the same colour, and by
# permuting colours it is enough to ask for both at colour 0.
# forced_same(i, j): adding the edge (i, j) leaves the graph uncolourable,
# asked through a selector so one solver answers every pair.
SEL = n * KC + 1
sel_pool = {}


def forced_diff(u, v):
    return not s.solve(assumptions=[1 + u * KC, 1 + v * KC])


def forced_same_pair(u, v):
    sv = sel_pool.get((u, v))
    if sv is None:
        sv = SEL + 1 + len(sel_pool)
        for c in range(KC):
            s.add_clause([-sv, -(1 + u * KC + c), -(1 + v * KC + c)])
        sel_pool[(u, v)] = sv
    return not s.solve(assumptions=[sv])


def ask(pairs, want_same):
    hit = 0
    for u, v in pairs[:SAMP]:
        hit += bool(forced_same_pair(u, v) if want_same else forced_diff(u, v))
    return hit


hs = ask(same, True)
print(f"  same-class pairs forced to agree:  {hs}/{SAMP}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
hd = ask(diff, False)
print(f"  cross-class pairs forced to differ: {hd}/{SAMP}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nrigidity: {100*(hs+hd)/(2*SAMP):.1f}% "
      f"(100% would mean uniquely {KC}-colourable, hence chi_c = {KC})",
      flush=True)
print("DONE", flush=True)
