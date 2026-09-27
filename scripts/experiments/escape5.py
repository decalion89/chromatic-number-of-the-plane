"""What would have to be forced for mu_5 to reach 3?  Ask the same instrument.

At four colours the escape analysis answered this completely: enumerate the
colourings that keep mu low, read off their monochromatic pairs inside N(p),
and compute the minimal family of pairs that meets them all.  That named the
requirement, and then showed it was the wrong SHAPE -- a disjunction of width
3-4 where the rotation stack needs width 2.

Nobody has asked it at five.  mu_5 = 2 means escapes exist that put only TWO
colours on N(p), and killing every one of them is exactly what raising mu_5 to
3 requires.  Each escape splits N(p) into two classes, so it carries many
monochromatic pairs -- the question is whether a small family of pairs meets
them all, and at what distances those pairs sit.

Three outcomes, all informative:

  * a small hitting set at reachable distances -- the first concrete recipe for
    moving mu at five, and the thing to build next;
  * a large one -- the requirement is real but diffuse, and no single gadget
    will do it;
  * an escape with no monochromatic non-adjacent pair -- then no family of
    forced pairs can raise mu_5 at that point at all, which is a genuine
    impossibility rather than a failure to find something.

Run on the small 5-chromatic graphs, where a colouring costs under a second.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
K = 5
for name in ("five_247_c.json", "five_tuned_1_1.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    print(f"\n  {name}  n={n} deg="
          f"{2.0*sum(len(a) for a in g.adj)/2/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
    rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
    S = set(g.vertices)
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    best = None
    for c in deg_order[:40]:
        for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
            for p in g.vertices:
                z = rot(p)
                if z in S:
                    continue
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
                if best is None or len(nb) > len(best[1]):
                    best = (z, sorted(nb))
    P0, NB = best
    print(f"    target |N(p)| = {len(NB)}   [{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    assert s.solve()
    # escapes at five: N(p) shows at most TWO colours, i.e. three are absent,
    # and colour symmetry lets those be 0, 1, 2
    escapes = []
    for _ in range(150):
        s.conf_budget(4_000_000)
        ass = [-X(u, c) for u in NB for c in (0, 1, 2)]
        if s.solve_limited(assumptions=ass) is not True:
            break
        pos = set(l for l in s.get_model() if l > 0)
        col = {u: next(c for c in range(K) if X(u, c) in pos) for u in NB}
        escapes.append(col)
        s.add_clause([-X(u, col[u]) for u in NB])
    s.delete()
    print(f"    {len(escapes)} escapes (two colours on N(p)) sampled"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if not escapes:
        print("    none -- mu_5 is already above 2 here!", flush=True)
        continue
    NBset = set(NB)
    adj = {(a, b) for a in NB for b in g.adj[a] if b in NBset}
    sets, freq = [], Counter()
    for col in escapes:
        mono = set()
        for i, a in enumerate(NB):
            for b in NB[i + 1:]:
                if col[a] == col[b] and (a, b) not in adj and (b, a) not in adj:
                    mono.add((a, b)); freq[(a, b)] += 1
        sets.append(mono)
    none_killable = sum(1 for m in sets if not m)
    print(f"    {none_killable} escapes have no monochromatic non-adjacent "
          f"pair; {len(freq)} distinct pairs appear", flush=True)
    if none_killable:
        print("    -> no family of forced pairs can raise mu_5 at this point",
              flush=True)
        continue
    unc = list(range(len(sets)))
    cover = []
    while unc:
        cnt = Counter()
        for i in unc:
            for pr in sets[i]:
                cnt[pr] += 1
        pr, _ = cnt.most_common(1)[0]
        cover.append(pr)
        unc = [i for i in unc if pr not in sets[i]]
    print(f"    a hitting set of {len(cover)} pairs meets every escape:",
          flush=True)
    for pr in cover[:10]:
        d2 = (g.vertices[pr[0]] - g.vertices[pr[1]]).norm2()
        print(f"       {pr}  d^2 = {d2}   ({float(d2):.4f})", flush=True)
    json.dump({"graph": name, "neighbourhood": NB, "escapes": len(escapes),
               "cover": [[a, b, str((g.vertices[a] - g.vertices[b]).norm2())]
                         for a, b in cover]},
              open(f"{ROOT}/data/escape5_{name}", "w"))
    print(f"    written data/escape5_{name}   [{time.time()-t0:.0f}s]",
          flush=True)
