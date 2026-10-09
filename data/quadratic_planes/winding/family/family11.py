"""family11.py d [k] [n] [--milp TL]: the configuration for d = 11 (mod 24).

U = G u G*u_1 u G*conj(u_1) u G*u_n, where G = rational unit vectors with denominator dividing 5^k,
u_1 = (1 + i sqrt d)^2/(1 + d) and u_n = (n + i sqrt d)^2/(n^2 + d), with n > 0, n = 3 (mod 4), n = 1 (mod 3^(s+1)),
s = v_3(d + 1); by default n = 1 + 2*3^(s+1).  notes/four_colours_11_mod_12.md (Theorem 1b) proves infeasibility when 5^k > (2 sqrt10/3) n (n + d).
Writes in_f11_{d}_{L}.json (format of certify_w2.py / check_w.py); with --milp runs the floating-point solver."""
import sys, json, math
from fractions import Fraction as Fr
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # winding/: kapparel
from family23 import rational_units


def v3(x):
    s = 0
    while x % 3 == 0: s += 1; x //= 3
    return s


def default_n(d):
    return 1 + 2 * 3 ** (v3(d + 1) + 1)


def proven_k(d, n):
    k = 0
    while 9 * 25 ** k <= 40 * (n * (n + d)) ** 2: k += 1     # 5^k > (2 sqrt10/3) n (n+d)  <=>  9 * 25^k > 40 (n(n+d))^2
    return k


def config(d, k, n):
    assert d % 24 == 11 and n > 0 and n % 4 == 3 and (n - 1) % 3 ** (v3(d + 1) + 1) == 0
    types = []
    for nn, sg in ((1, 1), (1, -1), (n, 1)):
        t = nn * nn + d
        types.append((Fr(nn * nn - d, t), Fr(2 * nn * sg, t)))        # vector (X, sqrt d * Y) with X real, Y real
    vecs = []
    for g1, g2 in rational_units(5 ** k):
        vecs.append((g1, Fr(0), g2, Fr(0)))
        for X, Y in types:                                            # gamma * (X + i Y sqrt d)
            vecs.append((g1 * X, -g2 * Y, g2 * X, g1 * Y))
    L = 1
    for v in vecs:
        for x in v: L = L * x.denominator // math.gcd(L, x.denominator)
    seen, U = set(), []
    for v in vecs:
        w = tuple(int(x * L) for x in v)
        key = max(w, tuple(-x for x in w))
        if key not in seen:
            seen.add(key); U.append(list(key))
    for a, b, c, e in U:
        assert a * a + d * b * b + c * c + d * e * e == L * L and a * b + c * e == 0
    return L, U


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    d = int(args[0]); n = int(args[2]) if len(args) > 2 else default_n(d)
    k = int(args[1]) if len(args) > 1 and args[1] != "-" else proven_k(d, n)
    L, U = config(d, k, n)
    path = f"in_f11_{d}_{L}.json"
    json.dump({"d": d, "D": L, "units": U}, open(path, "w"))
    print(f"d={d} n={n} k={k} N={5**k}: L={L}, {len(U)} vectors up to sign -> {path}; proven bound needs k >= {proven_k(d, n)}",
          flush=True)
    if "--milp" in sys.argv:
        import time
        from kapparel import relation_basis, solve_rel
        t0 = time.time(); ker = relation_basis(U)
        name, _ = solve_rel(len(U), ker, Fr(1, 3), 600)
        print(f"   MILP: {name} [{time.time()-t0:.0f}s]", flush=True)
