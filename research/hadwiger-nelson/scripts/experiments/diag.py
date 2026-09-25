import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/fastblock.py").read()
     .split('print(f"necklace')[0])
from hn.homcol import lattice_basis, on_lattice, has_homomorphism
t1 = time.time()
print(f"necklace {len(uniq)} points, {len(ints)} directions, "
      f"largest entry {max(max(abs(x) for x in v) for v in ints)} digits "
      f"{len(str(max(max(abs(x) for x in v) for v in ints)))}", flush=True)
B = lattice_basis(ints)
print(f"  Hermite: rank {len(B)}, largest entry {len(str(max(max(abs(x) for x in b) for b in B)))} digits"
      f"  [{time.time()-t1:.1f}s]", flush=True)
red = on_lattice(ints)
print(f"  coordinates: largest {len(str(max(max(abs(x) for x in v) for v in red)))} digits"
      f"  [{time.time()-t1:.1f}s]", flush=True)
for n in (2, 3, 4, 5):
    t2 = time.time()
    small = sorted({tuple(x % n for x in v) for v in red})
    zero = any(not any(v) for v in small)
    b = True if zero else has_homomorphism(small, n)[0] is None
    print(f"  n = {n}: {'blocks' if b else 'escapes'} "
          f"({len(small)} residues){' [zero vector]' if zero else ''}  "
          f"[{time.time()-t2:.1f}s]", flush=True)
