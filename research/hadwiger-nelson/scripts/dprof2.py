"""The decay curve, calibrated against a graph that actually forces.

Diversity comes from randomising the solver's saved phases, not from assuming
vertices into colours.  Cadical ignores set_phases and returns the same model
every time -- one distinct colouring in eight tries -- so it samples nothing;
glucose honours them and gives eight out of eight.  Assumption-based diversity
works but biases the sample, and on a tight graph most assumptions are
unsatisfiable and every retry is expensive.  Phases are free and unbiased.

Y forces a monochromatic pair in every proper 4-colouring, so Y's curve has a
floor: the forced pair agrees in every colouring that exists.  Whatever ratio
Y decays at is the shape of forcing.  G alone and the control -- two copies a
thousand apart, no cross edge possible -- give the free end.  Independent
pairs would decay at 1/k exactly.
"""
import sys, time, random
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Y, build_graph
from hn.geometry import DEGREY_FIELD as F
from hn.homcol import closable_distance
from pysat.solvers import Solver

MARKS = (1, 2, 3, 5, 10, 20, 40)
t0 = time.time()


def curve_of(name, n, edges, k, cand, cap=40, stall=20):
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable  [{time.time()-t0:.0f}s]",
              flush=True)
        sv.delete()
        return
    rng = random.Random(2718)
    surv, out, since = cand, [], 0
    for s in range(cap):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + i)
                       for i in range(n * k)])
        sv.solve()
        m = sv.get_model()
        col = [next(c for c in range(k) if m[w * k + c] > 0)
               for w in range(n)]
        before = len(surv)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        since = 0 if len(surv) < before else since + 1
        if s + 1 in MARKS:
            out.append((s + 1, len(surv)))
        if not surv or since >= stall:
            break
    out.append((s + 1, len(surv)))
    sv.delete()
    chk = Solver(name="cd19", bootstrap_with=cls)
    hits = [(i, j) for i, j in surv
            if not chk.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    chk.delete()
    # Consecutive marks only: the marks are not evenly spaced, and
    # averaging a two-sample step with one-sample steps understates the rate.
    r = [b / a for (s1, a), (s2, b) in zip(out, out[1:])
         if s2 == s1 + 1 and a]
    ratio = sum(r) / max(1, len(r))
    print(f"  {name}: {n} pts, {len(edges)} edges, {k} colours, "
          f"{len(cand)} pairs, decay {out}, ratio {ratio:.3f} "
          f"(free would be {1/k:.3f}), {len(hits)} FORCED"
          f"  [{time.time()-t0:.0f}s]", flush=True)


def pairs_of(zf):
    out, n = [], len(zf)
    for i in range(n):
        ai, bi = zf[i]
        for j in range(i + 1, n):
            v = (ai - zf[j][0]) ** 2 + (bi - zf[j][1]) ** 2
            if v > 36.0:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
                continue
            out.append((i, j))
    return out


ypts = build_Y(F)
yz = [(float(p.x), float(p.y)) for p in ypts]
ye = [(min(a, b), max(a, b)) for a, b in build_graph(ypts).edges()]
yp = pairs_of(yz)
print(f"Y: {len(yz)} vertices, {len(ye)} edges, {len(yp)} pairs  "
      f"[{time.time()-t0:.0f}s]", flush=True)
curve_of("Y at four colours (it forces)", len(yz), ye, 4, yp)
curve_of("Y at five colours", len(yz), ye, 5, yp)

pts = build_G(F, as_graph=False)
gz = [(float(p.x), float(p.y)) for p in pts]
ge = [(min(a, b), max(a, b)) for a, b in build_graph(pts).edges()]
n = len(gz)
gp = pairs_of(gz)
print(f"G: {n} points, {len(ge)} edges, {len(gp)} pairs  "
      f"[{time.time()-t0:.0f}s]", flush=True)
curve_of("G at five colours", n, ge, 5, gp)
curve_of("G at four colours", n, ge, 4, gp)

cz = gz + [(x + 1000.0, y2) for x, y2 in gz]
ce = ge + [(a + n, b + n) for a, b in ge]
curve_of("control: two copies 1000 apart", 2 * n, ce, 5, pairs_of(cz))
