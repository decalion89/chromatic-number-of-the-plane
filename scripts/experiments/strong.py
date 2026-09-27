"""Run the stronger gate above 5 on the sets that block at 5.

Blocking says nothing about the quotients Z/n for n > 5, where the Cayley
graph stops being complete and can well be 5-colourable.  If any modulus
yields a phi whose Cayley graph needs five colours or fewer, the graph is
5-colourable outright and blocking at 5 bought nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import denominator_29_directions, periodic_screen

d = denominator_29_directions()
print(f"{len(d)} denominator-29 directions, which block at n = 5", flush=True)
t0 = time.time()
for n in range(6, 21):
    r = periodic_screen(d, n, cap=5, rounds=60)
    if r.get("colourable"):
        print(f"  n = {n:2d}: *** 5-COLOURABLE *** Cayley chi = "
              f"{r['cayley_chi']} after {r['phis_tried']} phi  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        break
    tag = ("no homomorphism at all" if r["exhausted"] and not r["phis_tried"]
           else f"{r['phis_tried']} phi tried, all needing more than 5"
                + (", exhausted" if r["exhausted"] else ""))
    print(f"  n = {n:2d}: {tag}  [{time.time()-t0:.0f}s]", flush=True)
