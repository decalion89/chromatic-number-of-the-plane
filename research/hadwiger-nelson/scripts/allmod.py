"""Block at every modulus up to five, not just at five.

A coset colouring mod n gives chi <= n, so a 6-chromatic graph needs NO coset
colouring at n = 2, 3, 4 or 5.  The necklace blocks at 2, 3 and 5 and fails at
4 -- and the n = 8 "failure" found earlier was exactly that, since {0,4} is a
subgroup of Z/8 and avoiding it is a homomorphism to Z/4.  So the screen
sharpens and simplifies at the same time: four cheap tests, no Cayley
chromatic numbers needed.

The question for growth is whether the field has the material: do all 450
directions of its 75 zeta_6-orbits block at 4 as well?
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
exec(open("/tmp/claude-0/-home-user-darwin-50/"
          "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/orbblock.py")
     .read().split("phi, _ = has_homomorphism(iv, 5)")[0])
from hn.homcol import has_homomorphism

t0 = time.time()
print(f"{len(iv)} directions from {len(set(orbit.values()))} orbits", flush=True)
for n in (2, 3, 4, 5):
    phi, _ = has_homomorphism(iv, n)
    print(f"  n = {n}: " + ("BLOCKS" if phi is None else "coset colouring")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
