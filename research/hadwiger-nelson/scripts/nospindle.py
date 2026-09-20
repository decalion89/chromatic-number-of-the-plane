"""No Eisenstein spindle exists over Q(zeta_21), and that is a theorem.

A spindle needs two points forced to the SAME colour at k colours, at some
distance d, and a rotation rho with |1 - rho| = 1/d: then the rotated copy
puts them one apart while still equal, which is the contradiction.

In the Eisenstein lattice at three colours the forcing is exact, since the
lattice is uniquely 3-colourable: p and q take the same colour precisely when
p - q lies in the index-3 sublattice, so

    d^2 = N(p - q)  is divisible by 3,  say d^2 = 3m.

Then |1 - rho|^2 = 2 - 2 Re(rho) = 1/(3m) gives Re(rho) = (6m - 1)/(6m), and

    rho = ( (6m - 1) +- sqrt( -(12m - 1) ) ) / (6m).

So rho lies in Q(zeta_21) exactly when sqrt(-(12m-1)) does -- and the only
imaginary quadratic subfields there are Q(sqrt-3) and Q(sqrt-7). The condition
is therefore 12m - 1 = 3t^2 or 7t^2.

Modulo 12 that is impossible: 12m - 1 = 11, and 3t^2 is 0 or 3, while 7t^2
needs t^2 = 5 mod 12, which is not a square there -- squares mod 12 are
0, 1, 4 and 9.

m = 1 is the classical case: 12 - 1 = 11, the Moser rotation's sqrt(-11),
which is exactly the subfield Q(zeta_21) does not have.
"""
import sys
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from math import isqrt

print("squares mod 12:", sorted({t * t % 12 for t in range(12)}))
print("12m - 1 mod 12 :", (12 * 7 - 1) % 12)
bad3 = sorted({3 * t * t % 12 for t in range(12)})
bad7 = sorted({7 * t * t % 12 for t in range(12)})
print(f"  3t^2 mod 12 takes {bad3}; 7t^2 mod 12 takes {bad7}")
print(f"  11 in either? {11 in bad3 or 11 in bad7}")

print("\nchecked directly for the first thousand m:")
hits = []
for m in range(1, 1001):
    v = 12 * m - 1
    for d in (3, 7):
        if v % d:
            continue
        t = v // d
        if isqrt(t) ** 2 == t:
            hits.append((m, d, isqrt(t)))
print(f"  m with 12m-1 = 3t^2 or 7t^2: {hits if hits else 'NONE'}")
print(f"\n  m = 1 gives 11, the Moser rotation's sqrt(-11) -- and Q(zeta_21)")
print("  has only sqrt(-3) and sqrt(-7) among its imaginary quadratics.")
print("\nSo the lattice-doubling route to chi = 4 is closed over Q(zeta_21):")
print("no rotation of that field puts two same-coloured Eisenstein points")
print("one apart, which is why every union measured there stayed 3-chromatic.")
