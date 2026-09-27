"""Change the field, which changes what operations exist at all.

Every negative here is about one field.  De Grey chose Q(sqrt3, sqrt5, sqrt7,
sqrt11) for HIS problem, and the field decides which rings are closable, which
rotations exist, and therefore what the catalogue of operations even contains.
A richer field is not a bigger search of the same space; it is a different
space.

Two things change and both are measurable.  The DIRECTIONS: a field with more
radicals has more unit vectors, so the closure of a seed is denser and its
rings are more populated.  And the ROTATIONS: closable_distance asks for
sqrt(4D - 1) in the field, so every new square class opens rings that were
simply unavailable.

Built the same way, so the comparison is clean: the same seed, the same
dihedral closure, the same census question.  The seed is Sa's own S -- it lies
in every field that contains de Grey's -- and what changes is only which
points the closure reaches and which rings can then be bitten.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.field import Field
from hn.geometry import Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = 4
t0 = time.time()
BASE = (3, 5, 7, 11)
EXTRA = [(), (2,), (13,), (2, 13), (17,), (19,)]
for add in EXTRA:
    gens = tuple(sorted(set(BASE) | set(add)))
    try:
        K = Field(gens)
    except Exception as e:
        print(f"field {gens}: cannot build -- {e}", flush=True)
        continue
    half = K.rational(Fr(1, 2))
    rot60 = Rotation(half, K.sqrt(3) * half)
    # how many rings does this field make closable, of the first hundred?
    def closable(D):
        v = K.rational(4) * K.rational(D) - K.rational(1)
        c = v.c
        if any(x for x in c[1:]):
            return False
        q = c[0]
        if q <= 0:
            return False
        for r in K._prod:
            t = q / r
            n2, d2 = t.numerator, t.denominator
            rn, rd = int(round(n2 ** .5)), int(round(d2 ** .5))
            if rn * rn == n2 and rd * rd == d2:
                return True
        return False
    Ds = [Fr(x, y) for y in range(1, 10) for x in range(1, 10 * y + 1)]
    Ds = sorted({D for D in Ds if D != 1})
    nclos = sum(1 for D in Ds if closable(D))
    print(f"field {gens}: dimension {K.dim}, {nclos} of {len(Ds)} rings "
          f"closable  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
