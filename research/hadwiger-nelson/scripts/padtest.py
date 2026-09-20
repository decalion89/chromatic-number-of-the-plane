"""Does the greedy core builder find cores that padding PROVES exist?

On Sa u rot(Sa) at four colours rho = 5, with the forcing set
S = [107, 208, 502, 568, 618].  Padding then says: for any pivot p, the set
N(p) u (S minus N(p) minus p) contains S, so it forces, so T = S minus N(p)
minus p is a CORE of p of size at most 5.

That is a proof, not a search.  So if cegar_core cannot find a core at those
pivots, the greedy construction is unreliable -- and its verdict of "no core
in 60 steps" on the folded union at five colours says nothing at all.

Both are run here on the same pivots: what padding guarantees, and what the
search returns.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, cegar_core, is_core, pressure

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))

pts = list(degrey.build_Sa())
seen = set(pts)
for p in list(pts):
    q = ROT(p)
    if q not in seen:
        seen.add(q)
        pts.append(q)
g = build_graph(pts)
S = [107, 208, 502, 568, 618]
rel = ColourRelations(g, 4)
deg = g.degrees()
print(f"Sa u rot(Sa): {g.n} vertices; rho = 5 with S = {S}", flush=True)

for p in sorted(range(g.n), key=lambda v: -deg[v])[:4] + S[:2]:
    T = [v for v in S if v != p and v not in g.adj[p]]
    t = time.time()
    padded = is_core(rel, p, T)
    pr = pressure(rel, p)
    found, ok = cegar_core(rel, p, limit=40)
    print(f"  pivot {p:4} (deg {deg[p]:3}, pressure {pr}): "
          f"padding gives a core of {len(T)} -> verified {padded}; "
          f"cegar {'found ' + str(len(found)) if ok else 'FAILED in 40 steps'}"
          f"  [{time.time()-t:.0f}s]", flush=True)
