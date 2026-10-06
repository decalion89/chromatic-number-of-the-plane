"""shrink_core.py W.json CORE.cnf OUT.json: keep the vertices that occur in the clauses of a drat-trim core of the
certify4.py formula (variables x(v, c) = 4v + c + 1), and the listed cycles that lie inside them."""
import json, sys
W = json.load(open(sys.argv[1]))
keep = set()
for line in open(sys.argv[2]):
    if line.startswith(("p", "c")):
        continue
    for t in line.split():
        l = abs(int(t))
        if l:
            keep.add((l - 1) // 4)
pts = W["points"]; zero = pts.index([0] * len(pts[0]))
keep.add(zero)
order = sorted(keep); new = {v: i for i, v in enumerate(order)}
cyc = [[new[v] for v in c] for c in W["cycles"] if all(v in new for v in c)]
json.dump({"D": W["D"], "units": W["units"], "points": [pts[v] for v in order], "cycles": cyc}, open(sys.argv[3], "w"))
print(f"{len(pts)} -> {len(order)} points; {len(W['cycles'])} -> {len(cyc)} cycles")
