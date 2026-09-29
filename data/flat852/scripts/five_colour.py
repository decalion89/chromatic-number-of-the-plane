"""five_colour.py NAME.json -- find and verify a proper 5-colouring of the graph in NAME.json (tabucol, k = 5);
writes NAME.5colouring.json. Together with non-4-colourability this gives chromatic number exactly 5."""
import sys, json
import numpy as np
from tabu import tabucol

d = json.load(open(sys.argv[1]))
n = d["n"]; E = np.array([[a, b] for a, b, j in d["edges"]], dtype=np.int64)
for seed in range(1, 6):
    col = tabucol(n, E, k=5, maxiter=5_000_000, seed=seed)
    if col is not None:
        break
ok = col is not None and len(col) == n and all(0 <= c < 5 for c in col) and bool(np.all(col[E[:, 0]] != col[E[:, 1]]))
print(f"{sys.argv[1]}: proper 5-colouring found and verified on all {len(E)} edges: {ok}; colour class sizes "
      f"{np.bincount(col, minlength=5).tolist() if col is not None else None}")
if ok:
    json.dump({"graph": sys.argv[1], "colours": [int(c) for c in col]},
              open(sys.argv[1].replace(".json", ".5colouring.json"), "w"))
