import sys, json, time
import numpy as np
from flat import *
from lns import lns_repair
D, CD, U = directions()
P = np.load("g1023c_pts.npy"); st = json.load(open("g1023c_state.json"))
tri = st["triangle"]; col = st["col_partial"]
E, J = build_edges(P, U); n = len(P)
adj = [set() for _ in range(n)]
for a, b in E: adj[a].add(int(b)); adj[b].add(int(a))
new = [v for v in range(n) if col[v] < 0]
print(f"n={n}, new={len(new)}", flush=True)
t = time.time()
out = lns_repair(n, adj, col, new, tri, rmax=int(sys.argv[1]) if len(sys.argv) > 1 else 6, tl=900)
if out is not None:
    ok = all(out[a] != out[b] for a, b in E)
    print(f"LNS found a colouring (verified {ok}) in {time.time() - t:.0f}s", flush=True)
    json.dump({"triangle": tri, "col": out}, open("g1023c_colour.json", "w"))
else:
    print(f"LNS failed in {time.time() - t:.0f}s", flush=True)
