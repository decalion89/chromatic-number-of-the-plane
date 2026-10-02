"""mkcnf_state.py STATE OUT.cnf FORCE(a,b,c,e) [K]: the forcing formula of growforce4.py for a saved state:
proper K-colouring of the unit-distance graph on the state's points, plus c(0) != c(m), with that pair fixed to
colours 0, 1 (as kissat3 in growforce4.py). Writes OUT.cnf and OUT.meta.json (n, i0, im, edge count)."""
import json, sys
import numpy as np
st = json.load(open(sys.argv[1]))
out = sys.argv[2]
m = tuple(int(x) for x in sys.argv[3].split(","))
K = int(sys.argv[4]) if len(sys.argv) > 4 else 4
P = [tuple(p) for p in st["P"]]
U = [tuple(u) for u in st["U"]]
idx = {p: i for i, p in enumerate(P)}
assert len(idx) == len(P)
i0, im = idx[(0, 0, 0, 0)], idx[m]
E = [(i0, im)]
for i, p in enumerate(P):
    for u in U:
        j = idx.get((p[0] + u[0], p[1] + u[1], p[2] + u[2], p[3] + u[3]))
        if j is not None and i < j:
            E.append((i, j))
n = len(P)
with open(out, "w") as f:
    f.write(f"p cnf {K * n} {n + K * len(E) + 2}\n")
    f.write("".join(" ".join(str(K * v + c + 1) for c in range(K)) + " 0\n" for v in range(n)))
    f.write("".join(f"{-(K * a + c + 1)} {-(K * b + c + 1)} 0\n" for a, b in E for c in range(K)))
    f.write(f"{K * i0 + 1} 0\n{K * im + 2} 0\n")
json.dump({"n": n, "i0": i0, "im": im, "edges": len(E) - 1, "K": K, "m": m, "D": st["D"]}, open(out[:-4] + ".meta.json", "w"))
print(f"n={n} edges={len(E) - 1} i0={i0} im={im} D={st['D']} units={len(U)}")
