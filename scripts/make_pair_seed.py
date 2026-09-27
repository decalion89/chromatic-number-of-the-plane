"""Make a MODE=apart / MODE=same growth seed from a graph JSON: the same points, units and colouring,
with the target pair A, B set (grow_lean.py reads the keys A and B).
usage: make_pair_seed.py graph.json A B out.json"""
import sys, json
from fractions import Fraction as Fr
src, A, B, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
d = json.load(open(src))
assert 0 <= A < len(d["points"]) and 0 <= B < len(d["points"]) and A != B
d["A"], d["B"], d["status"] = A, B, "seed"
d["note"] = f"{src} with target pair A={A}, B={B}"
json.dump(d, open(out, "w"))
print(out, len(d["points"]), "points, pair", A, B)
