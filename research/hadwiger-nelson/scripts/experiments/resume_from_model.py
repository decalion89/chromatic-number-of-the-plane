"""Turn a solver model (kissat 'v' lines, or a DIMACS-style model file) for a saved pre_cdcl
instance into a checkpoint with a verified colouring, so the growth can resume from it."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
inst, model_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
d = json.load(open(inst)); K = 5
lits = []
for line in open(model_file):
    if line.startswith("v"): lits += [int(x) for x in line.split()[1:]]
    elif line.startswith("s"): print(line.strip())
pos = {l for l in lits if l > 0}
n = len(d["points"])
col = []
for v in range(n):
    cs = [c for c in range(K) if (1 + v * K + c) in pos]
    col.append(cs[0] if cs else 0)
# verify against the stored-unit edges
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]; U = [mk(xy) for xy in d["units"]]
key = lambda x, y: (round(x, 9), round(y, 9))
vk = {key(p.fx, p.fy): i for i, p in enumerate(V)}
bad = sum(1 for i, p in enumerate(V) for u in U for j in [vk.get(key(p.fx + u.fx, p.fy + u.fy))] if j is not None and col[i] == col[j])
print(f"colouring of {n} points: {bad} monochromatic unit steps")
assert bad == 0
d["colouring"] = col; d["status"] = "checkpoint"
json.dump(d, open(out, "w"))
print("written", out)
