"""Rebuild Z, record it exactly, and shrink it to a 4-critical core."""
import sys, time, json
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from hn.coloring import find_uncolorable_core, is_k_colorable

K1 = Field((3, 11, 247))
pts = build_Sa(K1)
rot = rotation_joining(Fr(1), K1).about(pts[25])
seen, H = set(), []
for p in pts:
    for q in (p, rot(p)):
        if q not in seen: seen.add(q); H.append(q)
rot2 = rotation_joining(Fr(64, 9), K1).about(H[157])
seen2, Z = set(), []
for p in H:
    for q in (p, rot2(p)):
        if q not in seen2: seen2.add(q); Z.append(q)
g = build_graph(Z)
print(f"Z: n={g.n}, m={sum(len(a) for a in g.adj)//2}", flush=True)

def dump(pt):
    return [[[c.numerator, c.denominator] for c in pt.x.c],
            [[c.numerator, c.denominator] for c in pt.y.c]]

json.dump({
    "field_generators": list(K1.gens),
    "basis": "coefficient i multiplies sqrt(product of generators selected by bit i)",
    "recipe": [
        "S = de Grey's 39 points (hn.degrey.S_POINTS), taken in Q(sqrt3,sqrt11,sqrt247)",
        "Sa = orbit of S under the 12-element group <rot60, reflection in x>",
        "H  = Sa u rho1(Sa), rho1 = 60-degree rotation about Sa[25] (a vertex, not the origin)",
        "     glue circle: the 20 points at distance 1 from Sa[25]; it is NOT capped",
        "     H has 8 forced-equal pairs at 4 colours; Sa itself has none",
        "Z  = H u rho2(H), rho2 = rotation about H[157] by 2*arcsin(3/16),",
        "     cos 119/128, sin 384*sqrt(247)/16384, the spindle at |H[157]-H[327]| = 8/3",
    ],
    "n": g.n,
    "points": [dump(p) for p in g.vertices],
}, open("/home/user/darwin-50/research/hadwiger-nelson/data/five_247.json", "w"))
print("  written data/five_247.json", flush=True)

t0 = time.time()
core = find_uncolorable_core(g, 4, timeout=1200, rounds=80)
print(f"\n4-critical-ish core: n={core.n}, "
      f"m={sum(len(a) for a in core.adj)//2}   [{time.time()-t0:.0f}s]", flush=True)
ok, _ = is_k_colorable(core, 4, timeout=1200)
print(f"  core is 4-colourable: {ok}", flush=True)
json.dump({"field_generators": list(K1.gens), "n": core.n,
           "points": [dump(p) for p in core.vertices]},
          open("/home/user/darwin-50/research/hadwiger-nelson/data/five_247_core.json", "w"))
print("  written data/five_247_core.json", flush=True)
