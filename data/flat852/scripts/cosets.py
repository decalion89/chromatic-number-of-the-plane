"""cosets.py -- class of a point of M = L + conj(L) in M/O = F7*omega + F7*conj(omega) (O = Z[zeta21]).
A point x with integer vector c = 7x: c mod 7 = a*(7 omega mod 7) + b*(7 conj(omega) mod 7), (a, b) in F7^2.
x lies in the L-coset indexed by b (rows) and the conj(L)-coset indexed by a (columns)."""
import sys
import numpy as np
from flat import *

W1 = np.array(OM7) % 7
W2 = np.array(conj(tuple(OM7))) % 7
TAB = {}
for a in range(7):
    for b in range(7):
        TAB[tuple((a * W1 + b * W2) % 7)] = (a, b)
assert len(TAB) == 49


def cls(P):
    out = []
    for c in np.asarray(P) % 7:
        out.append(TAB.get(tuple(c), None))
    return out


if __name__ == "__main__":
    P = np.load(sys.argv[1])
    C = cls(P)
    assert all(c is not None for c in C), "point not in M"
    M = np.zeros((7, 7), dtype=int)
    for a, b in C:
        M[b, a] += 1
    print(f"{len(P)} points; rows b (L-cosets) x columns a (conj(L)-cosets):")
    for b in range(7):
        print("  b=%d: " % b + " ".join("%5d" % x for x in M[b]))
    D, CD, U = directions()
    # check: directions -- even: (0,0); omega*mu42: (a,0); conj(omega)*mu42: (0,b)
    cu = cls(U)
    print("direction classes: mu42", set(cu[j] for j in range(0, 84, 2)), " omega mu42 b-values",
          set(c[1] for c in (cu[j] for j in range(1, 84, 2))), " conj(omega) mu42 a-values", set(c[0] for c in cu[84:]))
