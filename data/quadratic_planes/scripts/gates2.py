"""For d = 11 (mod 12) and a denominator D: how many of the unit vectors with denominator dividing D are NOT
integral at each 'gate', i.e. a place of K(i) whose residue graph is 2- or 3-colourable:
  2-gate  (d = 7 mod 8): places above 2 with residue field F_2, residue graph a perfect matching (2 colours);
  3-gate  (always, as d = 2 mod 3): places above 3 with residue F_9, split in K(i)/K: H_9 (3 colours);
  5-gate  (d = 0, 1, 4 mod 5): places above 5 with residue F_5, split in K(i)/K: H_5 (3 colours).
If some gate has count 0, reduction there colours every graph built from these directions."""
import sys
from units_fast import units_fast
from gates import vp, padic_sqrt, val_via_minus_d, val_split_gauss, val_split_complete


def sqrt_mod_2k(a, k):
    """t with t^2 = a mod 2^k, a = 1 mod 8"""
    t = 1
    for j in range(3, k):
        if (t * t - a) % (2 ** (j + 1)):
            t += 2 ** (j - 1)
    return t % 2 ** k


def val2(a, b, c, e, D, d, sign=1):
    k = 40
    M = 2 ** k
    t = sign * sqrt_mod_2k((-d) % M, k)
    X = a * a + c * c - d * (b * b + e * e)
    Y = 2 * (a * e - b * c)
    N = (X + Y * t) % M
    return (vp(N, 2) if N else 99) - 2 * vp(D, 2)


def report(d, D):
    U = units_fast(d, D)
    res = {"n": len(U)}
    if d % 8 == 7:
        res["2"] = sum(1 for u in U if val2(*u, D, d) != 0)
    res["3"] = sum(1 for u in U if val_via_minus_d(*u, D, d, 3) != 0)
    if d % 5 == 0:
        res["5"] = sum(1 for u in U if val_split_gauss(*u, D, d, 5) != 0)
    elif d % 5 in (1, 4):
        res["5"] = sum(1 for u in U if val_split_complete(*u, D, d, 5) != 0)
    return res


if __name__ == "__main__":
    for pair in sys.argv[1:]:
        d, D = map(int, pair.split(":"))
        r = report(d, D)
        closed = [g for g in ("2", "3", "5") if r.get(g) == 0]
        print(f"d={d} D={D}: {r}  {'CLOSED ' + ','.join(closed) if closed else 'all open'}", flush=True)
