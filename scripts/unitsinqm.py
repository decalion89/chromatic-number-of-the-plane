import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, json
from fractions import Fraction as Fr
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
ROOT = HN_DIR
for name in sys.argv[1:]:
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    E = edge_vectors(build_graph(P))
    B = echelon(E); C = [coords(B, v) for v in E]; r = len(B)
    bad = []
    for q in range(2, 41):
        k = sum(1 for c in C if all(x % q == 0 for x in c))
        if k: bad.append((q, k))
    from math import gcd
    g = [abs(__import__('functools').reduce(gcd, c)) for c in C]
    print(f"{name}: {len(E)} directions, rank {r}; q with some unit in qM (q, #units): {bad}")
    print(f"   content (gcd of coordinates) of each unit, histogram: " + str(sorted(__import__('collections').Counter(g).items())))
