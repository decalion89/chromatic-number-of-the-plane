"""The requirement is per point.  Which point asks for the least?

Raising mu_5 above 2 at a point p costs a family of pairs forced apart -- the
ones that meet every escape.  At one point of the 803-graph that family is a
single pair at squared distance 3.  Nothing says every point is so cheap, or
that 3 is the only distance on offer, and the cheapest requirement in the
graph is the one worth building for.

So: sweep the candidate points, compute each requirement to exhaustion (assume
the cover so far, ask for another escape, extend, repeat until none survives),
and tabulate.  What comes out is a menu -- how many pairs, at what distances --
and the entry with one pair at the friendliest distance is the target.

Three things make it cheap.  Colour symmetry turns "two colours on N(p)" into
"colours 0, 1, 2 absent", one assumption list.  Escapes come from one warm
solver with a blocking clause each.  And the 803-point graph colours in under
a second, so the whole sweep is solver calls, not solver waits.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
K = 5
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
S = set(g.vertices)
cells = defaultdict(list)
for i, p in enumerate(g.vertices):
    cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
cand = set()
for c in deg_order[:40]:
    for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
        for p in g.vertices:
            z = rot(p)
            if z not in S:
                cand.add(z)
scored, seenn = [], set()
for z in cand:
    zx, zy = z.fx, z.fy
    cx, cy = int(zx // 1), int(zy // 1)
    nb = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for i in cells.get((cx + dx, cy + dy), ()):
                ex, ey = gx[i] - zx, gy[i] - zy
                if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                   (g.vertices[i] - z).norm2() == one:
                    nb.append(i)
    if len(nb) >= 8:
        t = tuple(sorted(nb))
        if t not in seenn:
            seenn.add(t); scored.append((len(nb), t))
scored.sort(key=lambda u: -u[0])
print(f"  {len(scored)} neighbourhoods >= 8, largest {scored[0][0]}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            base.append([-X(v, a), -X(v, b)])
for x, y in g.edges():
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

menu = []
for sz, NB in scored[:30]:
    NB = list(NB); NBset = set(NB)
    adj = {(a, b) for a in NB for b in g.adj[a] if b in NBset}
    cover, ok = [], False
    for rounds in range(12):
        cnf = list(base)
        for (a, b) in cover:
            for c in range(K):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        if not s.solve():
            s.delete(); cover = None; break
        got = []
        for _ in range(60):
            s.conf_budget(2_000_000)
            if s.solve_limited(assumptions=[-X(u, c) for u in NB
                                            for c in (0, 1, 2)]) is not True:
                break
            pos = set(l for l in s.get_model() if l > 0)
            got.append({u: next(c for c in range(K) if X(u, c) in pos)
                        for u in NB})
            s.add_clause([-X(u, got[-1][u]) for u in NB])
        s.delete()
        if not got:
            ok = True; break
        sets = [{(a, b) for i, a in enumerate(NB) for b in NB[i+1:]
                 if col[a] == col[b] and (a, b) not in adj and (b, a) not in adj}
                for col in got]
        if any(not m for m in sets):
            cover = None; break
        unc = list(range(len(sets)))
        while unc:
            cnt = Counter()
            for i in unc:
                for pr in sets[i]:
                    cnt[pr] += 1
            pr, _ = cnt.most_common(1)[0]
            cover.append(pr)
            unc = [i for i in unc if pr not in sets[i]]
    if cover is None or not ok:
        print(f"    |N|={sz}: no finite pair requirement", flush=True)
        continue
    ds = [str((g.vertices[a] - g.vertices[b]).norm2()) for a, b in cover]
    menu.append((len(cover), sz, cover, ds))
    print(f"    |N|={sz}: {len(cover)} pair(s), d^2 = {ds}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
menu.sort()
print("\n  cheapest requirements first:", flush=True)
for k, sz, cov, ds in menu[:10]:
    print(f"    {k} pair(s)  |N|={sz}  {list(zip(cov, ds))}", flush=True)
dist = Counter(x for _, _, _, ds in menu for x in ds)
print(f"\n  distances asked for: {dict(dist)}", flush=True)
json.dump([{"pairs": k, "N": sz, "cover": [[a, b] for a, b in cov],
            "d2": ds} for k, sz, cov, ds in menu],
          open(f"{ROOT}/data/mu5_menu.json", "w"))
print("  written data/mu5_menu.json", flush=True)
