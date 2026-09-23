import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
exec(open("/home/user/darwin-50/research/hadwiger-nelson/scripts/gate.py").read().split("def gate(g, label):")[0])
from hn.homcol import periodic_screen, cayley_chromatic
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
NAME = sys.argv[1]; MODS = [int(x) for x in sys.argv[2:]]
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
E = edge_vectors(g); B = echelon(E); C = [coords(B, v) for v in E]
print(f"{NAME}: {len(E)} directions, module rank {len(B)}", flush=True)
t0 = time.time()
for n in MODS:
    r = periodic_screen(C, n, cap=5, rounds=300)
    print(f"  Z/{n}: {r}   [{time.time()-t0:.0f}s]", flush=True)
