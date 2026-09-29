"""task1_model.py -- exact checks of the 126 directions and of the mod-2 claim."""
import json, itertools
from fractions import Fraction as Fr
import numpy as np
from flat import *

D, CD, U = directions()
print(f"|D| = {len(D)}, |conj D| = {len(CD)}, |U| = {len(U)}")
assert len(set(D)) == 84 and len(set(CD)) == 84 and len(U) == 126
print(f"|D & conj D| = {len(set(D) & set(CD))} (should be the 42 of mu42)")

# omega exactly: minpoly 7x^12 - 13x^6 + 7
w = list(OM7)          # 7*omega
p = [49 * 0 + (1 if i == 0 else 0) * 1 for i in range(N)]
# compute omega^6 and omega^12 as Fractions via integer arithmetic: (7w)^6 = 7^6 w^6
w6 = [1] + [0] * 11
for _ in range(6):
    w6 = mul(w6, w)
w12 = mul(w6, w6)
# 7*w^12 - 13*w^6 + 7 = 0  <=>  7*(7w)^12 - 13*7^6*(7w)^6 + 7*7^12 = 0
lhs = [7 * a - 13 * 7 ** 6 * b + (7 * 7 ** 12 if i == 0 else 0) for i, (a, b) in enumerate(zip(w12, w6))]
print("7 w^12 - 13 w^6 + 7 == 0 exactly:", all(x == 0 for x in lhs))

# numerics of omega vs the formula zeta3/(zeta7 - zeta7^-1) - zeta3^2/(zeta7^2 - zeta7^-2)
z3, z7 = cmath.exp(2j * math.pi / 3), cmath.exp(2j * math.pi / 7)
om_formula = z3 / (z7 - 1 / z7) - z3 ** 2 / (z7 ** 2 - z7 ** -2)
print("|omega - formula| =", abs(cval(OM7) - om_formula))

# u * conj(u) = 1 exactly, i.e. (7u)(7 conj u) = 49 mod Phi21
ok_unit = all(mul(list(u), list(conj(u))) == [49] + [0] * 11 for u in U)
print("u*conj(u) == 1 exactly for all 126:", ok_unit)
# closed under negation
S = set(U)
print("closed under negation:", all(tuple(-x for x in u) in S for u in U))
# conj(U) == U
print("conj(U) == U:", set(conj(u) for u in U) == S)
# distinct as complex numbers, all of modulus 1; angles
vals = [cval(u) for u in U]
print("max ||u|-1| (float):", max(abs(abs(v) - 1) for v in vals))
md = min(abs(a - b) for a, b in itertools.combinations(vals, 2))
print("min distance between two directions (float):", md)
# all coefficients integral after scaling by 7; which directions lie in Z[zeta21] (all coords divisible by 7)
inO = [u for u in U if all(c % 7 == 0 for c in u)]
print("directions in Z[zeta21]:", len(inO), "(mu42)")

# ---------------- mod 2 ----------------
def res(u):
    return sum(((c % 2) << k) for k, c in enumerate(u))    # 1/7 read as 1: numerator mod 2
par = lambda x: bin(x).count("1") & 1

def solve_affine(rows, rhs):
    """number of b in F2^12 with par(b & r) = rhs for all rows"""
    piv = {}
    for r, t in zip(rows, rhs):
        for pb in sorted(piv, reverse=True):
            if (r >> pb) & 1:
                r ^= piv[pb][0]; t ^= piv[pb][1]
        if r == 0:
            if t:
                return 0
            continue
        pb = r.bit_length() - 1
        # keep rows reduced
        for q in list(piv):
            if (piv[q][0] >> pb) & 1:
                piv[q] = (piv[q][0] ^ r, piv[q][1] ^ t)
        piv[pb] = (r, t)
    return 2 ** (12 - len(piv))

def count_pairs(dirs):
    R = sorted(set(res(u) for u in dirs))
    tot = 0
    for a in range(4096):
        Z = [r for r in R if par(a & r) == 0]
        tot += solve_affine(Z, [1] * len(Z))
    return R, tot

def brute_pairs(dirs):
    """independent brute force over all (a,b) using numpy (4096 x 4096 pairs)"""
    R = np.array(sorted(set(res(u) for u in dirs)), dtype=np.int64)
    bits = ((R[:, None] >> np.arange(12)) & 1).astype(np.uint8)          # (|R|,12)
    A = ((np.arange(4096)[:, None] >> np.arange(12)) & 1).astype(np.uint8)  # (4096,12)
    dots = (A.astype(np.int64) @ bits.T.astype(np.int64)) & 1             # (4096,|R|): a.r
    # pair (a,b) good iff for every r: dots[a,r] | dots[b,r]
    zero = (dots == 0)                                                    # a.r == 0
    # (a,b) bad iff exists r with zero[a,r] & zero[b,r]
    Zf = zero.astype(np.float32)
    bad = (Zf @ Zf.T) > 0
    return int((~bad).sum())

for name, dirs in (("D (84)", D), ("conj D (84)", CD), ("U (126)", U)):
    R, tot = count_pairs(dirs)
    bf = brute_pairs(dirs)
    print(f"{name}: {len(R)} distinct residues mod 2 (zero residue present: {0 in R}); "
          f"ordered pairs (a,b) with (a.r,b.r)!=(0,0) for all r: {tot} (brute force: {bf}); "
          f"row spaces {{a,b,a+b}}: {tot / 6:g}")

out = {"convention": "integer vector v means (v_0 + v_1 z + ... + v_11 z^11)/7, z = zeta21 = exp(2 pi i/21)",
       "D": [list(u) for u in D], "conjD": [list(u) for u in CD], "U": [list(u) for u in U]}
json.dump(out, open("directions.json", "w"))
print("wrote directions.json")
