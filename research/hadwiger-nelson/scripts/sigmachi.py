"""The LOCAL chromatic number after sigma -- one triangle can never help, two can.

A single ring triangle is capped at three, and the proof is an angle.  A point u
of N(h) reaches the ring only at the angles phi +- theta0, so it has at most two
ring neighbours, and those two are 2*theta0 = 146.44 degrees apart -- never a
multiple of 120, so never in the same coset, so never in the same triangle.  So
u has AT MOST ONE neighbour in any given triangle, keeps a list of two of its
three colours, and N(h) is a union of paths and even cycles, which are
2-choosable.  mu(N(h) u T) <= 3 always.  No graph, no enrichment, no search.

Two triangles is a different matter.  There u may be adjacent to one vertex of
each, and if those carry different colours its list drops to ONE -- and then the
internal edges of N(h) have no slack left.  That is the first configuration in
this project where the local floor can move at all, and sigma is what builds it:
one round on the 803-graph gives 8323 points with 86 896 cross-ring triangles,
and its best hub has 54 neighbours, a ring of 42, and 30 neighbours tied to the
ring twice over.

So compute chi of N(h) u ring(h) directly -- purely local, no ambient colouring
needed, a hundred vertices at most.  Before sigma it was 3 at every hub of every
graph here.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 1
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
r11 = F.sqrt(11); sixth = F.rational(Fr(1, 6)); third = F.rational(Fr(1, 3))
one = F.rational(1)
def sig(p, c, conj):
    s = -r11 if conj else r11
    dx, dy = p.x - c.x, p.y - c.y
    return Point(c.x + (dx - s * dy) * sixth, c.y + (s * dx + dy) * sixth)
pts = {}
for x, y in d["points"]:
    q = Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
    pts[(round(q.fx, 9), round(q.fy, 9))] = q
for _ in range(ROUNDS):
    G0 = build_graph(list(pts.values()))
    for c, u in G0.edges():
        for (a, b) in ((c, u), (u, c)):
            for conj in (False, True):
                p = sig(G0.vertices[b], G0.vertices[a], conj)
                k = (round(p.fx, 9), round(p.fy, 9))
                if k not in pts: pts[k] = p
G = build_graph(list(pts.values())); n = G.n
adj = defaultdict(set)
for x, y in G.edges():
    adj[x].add(y); adj[y].add(x)
hx = [q.fx for q in G.vertices]; hy = [q.fy for q in G.vertices]
print(f"sigma^{ROUNDS} of {NAME}: n={n} edges={sum(len(a) for a in adj.values())//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
ring = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if abs(dd - 1/3) < 1e-9 and (G.vertices[i]-G.vertices[j]).norm2() == third:
            ring[i].append(j); ring[j].append(i)
print(f"  {len(ring)} hubs, largest ring {max(len(v) for v in ring.values())}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def chi(verts, edges, cap=6):
    for k in range(1, cap + 1):
        idx = {v: i for i, v in enumerate(verts)}
        X = lambda v, c: 1 + idx[v] * k + c
        cnf = [[X(v, c) for c in range(k)] for v in verts]
        for a, b in edges:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if ok: return k
    return cap + 1

tally = Counter(); best = None
order = sorted(ring, key=lambda h: -(len(adj[h]) + len(ring[h])))
for h in order[:1500]:
    V = sorted(set(adj[h]) | set(ring[h]))
    if len(V) < 6: continue
    S = set(V)
    Eloc = [(a, b) for a in V for b in adj[a] if b in S and b > a]
    k = chi(V, Eloc)
    tally[k] += 1
    if best is None or k > best[0]:
        best = (k, h, len(adj[h]), len(ring[h]), len(V), len(Eloc))
        print(f"    chi = {k} at h={h}: |N|={len(adj[h])} |ring|={len(ring[h])} "
              f"|X|={len(V)} edges={len(Eloc)}   [{time.time()-t0:.0f}s]", flush=True)
        if k >= 4:
            json.dump({"source": NAME, "rounds": ROUNDS, "hub": h,
                       "N": sorted(adj[h]), "ring": ring[h], "chi": k},
                      open(f"{ROOT}/data/sigma_localchi.json", "w"))
print(f"\n  chi(N(h) u ring) over {sum(tally.values())} hubs: "
      f"{dict(sorted(tally.items()))}   [{time.time()-t0:.0f}s]", flush=True)
