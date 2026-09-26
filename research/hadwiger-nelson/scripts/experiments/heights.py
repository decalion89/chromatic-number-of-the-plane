"""How tall are the available unit steps?

The necklace's direction coordinates reach 87 digits, and that is what makes
the Hermite reduction -- hence the honest blocking test -- cost minutes.  The
height comes from the steps: each is a Hilbert-90 quotient a/abar and a tall a
gives a tall step.  If short ones are plentiful the necklace can be rebuilt
inside them and the whole question becomes cheap.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import gcd
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/hn/"
           "kneck2.py").read()
head = src[:src.index("pairsum = {}")]
exec(head)


def height(u):
    h = 1
    for part in u:
        for q in part:
            h = max(h, abs(q.numerator), q.denominator)
    return len(str(h))


hs = Counter(height(u) for u in steps)
print(f"\nstep heights in digits: "
      f"{sorted((d, c) for d, c in hs.items())}  [{time.time()-t0:.0f}s]",
      flush=True)
short = [u for u in steps if height(u) <= 4]
print(f"{len(short)} steps of at most 4 digits, "
      f"{len([u for u in steps if height(u) <= 2])} of at most 2",
      flush=True)
orb = set()
for u in short:
    v, best = u, None
    for _ in range(6):
        k = tuple(tuple(q for q in part) for part in v)
        if best is None or k < best:
            best = k
        v = kmul(v, Z6)
    orb.add(best)
print(f"  they fall into {len(orb)} zeta_6-orbits", flush=True)
print(f"  the 44 steps of the built necklace have heights "
      f"{sorted(Counter(height(u) for u in __import__('pickle').load(open('/tmp/hn/kneck.pkl','rb'))).items())}",
      flush=True)
