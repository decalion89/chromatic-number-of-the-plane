"""How far each graph is from uniquely k-colourable, in one number.

A uniquely k-colourable graph has pressure exactly k-1 at EVERY vertex: if one
had two free colours, switching between them would move it to a different
class and give a genuinely different partition. Pressure is at most k-1
always, so the defect k-1-pressure(v) is a per-vertex distance to the target,
and it is monotone in added points -- something a search can climb, unlike a
yes/no verdict.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_S, build_Sa, build_Y
from hn.forced import ColourRelations, unique_colouring_defect
from hn.graph import build_graph
from hn.mixed import three_hexagon_gadget

cases = [("S", build_S(), (3, 4)), ("Sa", build_Sa(), (4,)),
         ("Y", build_Y(), (4,)), ("G", build_G(), (5,)),
         ("3-hexagon gadget", three_hexagon_gadget()[2], (4, 5))]
t0 = time.time()
for name, obj, ks in cases:
    g = obj if hasattr(obj, "vertices") else build_graph(obj)
    for k in ks:
        rel = ColourRelations(g, k)
        if not rel.colourable:
            rel.close()
            continue
        d = unique_colouring_defect(rel)
        rel.close()
        print(f"{name:>18} n={g.n:5d} k={k}: defects "
              f"{d['defect_histogram']}, at k-1: {d['at_k_minus_one']}/"
              f"{d['vertices']}  [{time.time()-t0:.0f}s]", flush=True)
