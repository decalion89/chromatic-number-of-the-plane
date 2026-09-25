"""Put the 4-critical necklace through the stronger gate.

Blocking at n = 5 is the bottom rung.  The necklace's directions block there,
which is what makes it the first 4-critical blocked unit-distance graph -- but
the condition a 6-chromatic graph really has to meet is that EVERY finite
abelian quotient has Cayley chromatic number at least 6.  Above 5 the Cayley
graph stops being complete, so its chromatic number has to be computed, and a
single hit anywhere would prove the direction set 5-colourable and the whole
substrate useless.

Run at every modulus from 5 up, by the CEGAR loop: solve for a homomorphism,
colour its Cayley graph, exclude it if that needs more than five colours.
Exhausting with no homomorphism at all is blocking, for that modulus.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, pickle, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import periodic_screen

with open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/neckdirs.pkl",
          "rb") as fh:
    dirs = [tuple(d) for d in pickle.load(fh)]
print(f"{len(dirs)} necklace directions, dimension {len(dirs[0])}", flush=True)

t0 = time.time()
for n in range(5, 21):
    r = periodic_screen(dirs, n, cap=5, rounds=60)
    if r.get("colourable"):
        print(f"  n = {n:2d}: *** 5-COLOURABLE *** Cayley chi = "
              f"{r['cayley_chi']} after {r['phis_tried']} phi  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        break
    tag = ("no homomorphism at all -- blocking"
           if r["exhausted"] and not r["phis_tried"]
           else f"{r['phis_tried']} phi tried, all needing more than 5"
                + (", exhausted" if r["exhausted"] else ""))
    print(f"  n = {n:2d}: {tag}  [{time.time()-t0:.0f}s]", flush=True)
