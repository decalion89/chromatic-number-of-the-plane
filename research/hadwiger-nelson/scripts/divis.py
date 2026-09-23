"""How divisible is each unit direction inside the edge module?  e in gM iff g
divides every coordinate of e in a Z-basis of M."""
import sys, json
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
for NAME in sys.argv[1:]:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P)
    E = edge_vectors(g); B = echelon(E); C = [coords(B, v) for v in E]
    divs = sorted((reduce(gcd, [abs(t) for t in c]), k) for k, c in enumerate(C))
    from collections import Counter
    print(f"{NAME}: {len(E)} directions, rank {len(B)}; divisibility g (e in gM) counts: {sorted(Counter(g for g, k in divs).items())}")
