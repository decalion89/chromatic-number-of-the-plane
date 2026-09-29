"""core_vertices.py PTS.npy CORE.cnf TRI OUTNAME -- vertices whose at-least-one clause lies in a drat-trim clause core
(plus the fixed triangle), then the 4-core of the induced subgraph. The induced subgraph is again not 4-colourable:
a colouring of it, extended by 'no colour' on the other vertices, would satisfy the core.
Writes OUTNAME_pts.npy and OUTNAME_state.json (triangle indices in the new order)."""
import sys, json
import numpy as np
from flat import *

D, CD, U = directions()
P = np.load(sys.argv[1]); corep = sys.argv[2]; tri = [int(x) for x in sys.argv[3].split(",")]; out = sys.argv[4]
n = len(P)
keep = set(tri)
nalo = 0
for line in open(corep):
    if line.startswith(("p", "c")):
        continue
    lits = [int(x) for x in line.split()[:-1]]
    if len(lits) == 4 and all(l > 0 for l in lits):
        v = (lits[0] - 1) // 4
        assert sorted(lits) == [4 * v + 1 + c for c in range(4)]
        keep.add(v); nalo += 1
E, J = build_edges(P, U)
adj = [set() for _ in range(n)]
for a, b in E:
    adj[a].add(int(b)); adj[b].add(int(a))
S = set(keep)
ch = True
while ch:
    ch = False
    for v in list(S):
        if v not in tri and len(adj[v] & S) < 4:
            S.discard(v); ch = True
L = sorted(S)
np.save(f"{out}_pts.npy", P[L])
json.dump({"triangle": [L.index(v) for v in tri], "source": sys.argv[1]}, open(f"{out}_state.json", "w"))
print(f"core ALO clauses: {nalo}; vertices kept (with triangle): {len(keep)}; after 4-core: {len(L)}; "
      f"triangle in new order: {[L.index(v) for v in tri]}")
