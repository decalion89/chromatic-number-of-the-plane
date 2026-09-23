"""The ceiling at four is a RAINBOW, and the geometry of one is readable.

five_247_c is 5-chromatic and vertex-critical, so every vertex p is a ceiling
point of H = G - p: mu_4(H, p) = 4, exactly, with no solving.  Two of the 803
have degree 4.  There mu_4 = 4 = |N(p)|: the cap is attained, which means the
four neighbours take four DIFFERENT colours in every 4-colouring of H -- six
pairs, pairwise forced apart, inside one unit circle.

Which of the six are edges (free, at 60 degrees) and which are genuine
forcings?  That is the whole content of the configuration, and it is a table
of squared distances.
"""
import sys, json
from fractions import Fraction as Fr
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
adj = [set() for _ in range(n)]
for x, y in g.edges():
    adj[x].add(y); adj[y].add(x)

for lim in (4, 5, 6):
    small = [p for p in range(n) if len(adj[p]) == lim]
    if not small: continue
    print(f"\n=== degree {lim}: {len(small)} ceiling points ===", flush=True)
    for p in small[:4]:
        nb = sorted(adj[p])
        e = [(u, v) for u, v in combinations(nb, 2) if v in adj[u]]
        print(f"  p={p}  N={nb}  edges inside N: {len(e)} of {len(nb)*(len(nb)-1)//2}", flush=True)
        for u, v in combinations(nb, 2):
            d2 = (g.vertices[u] - g.vertices[v]).norm2()
            tag = "EDGE (60 deg, free)" if v in adj[u] else "forced apart by H"
            print(f"     |{u}-{v}|^2 = {float(d2):.6f}   {tag}", flush=True)
