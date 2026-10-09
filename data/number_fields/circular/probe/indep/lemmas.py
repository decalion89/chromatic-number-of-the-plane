"""Referee checks of Lemmas 0, 1, 2, C, Q and of the base case k = 1 (exact arithmetic)."""
from fractions import Fraction as Fr
from itertools import product
from hgeom import Poly, funcs, ev, fl, gmul, gconj, gadd, gsub, rho_pow
from lp3x import lp_min
from sn import Pk, Yk, QN, cstar

I = (Fr(0), Fr(1))
H = (Fr(1, 2), Fr(1, 2))


def gpow(z, n):
    out = (Fr(1), Fr(0))
    for _ in range(n):
        out = gmul(out, z)
    return out


def is_gint(z):
    return z[0].denominator == 1 and z[1].denominator == 1


def lattice_pm(sign, box=Fr(3, 2)):
    """points of Lambda_sign \\ Z[i] = ((2 + sign i)/5) Z[i] minus Z[i] with both coordinates in [-box, box]"""
    g = (Fr(2, 5), Fr(sign, 5))
    out = []
    R = int(4 * box) + 4
    for a in range(-R, R + 1):
        for b in range(-R, R + 1):
            z = gmul(g, (Fr(a), Fr(b)))
            if abs(z[0]) <= box and abs(z[1]) <= box and not is_gint(z):
                out.append(z)
    return sorted(set(out))


print("== Lemma 1 (values at type points), k = 1..9, all |j| <= k, all q")
ok = True
for k in range(1, 10):
    N = 5 ** k
    for j in range(-k, k + 1):
        rj = rho_pow(j)
        v = gmul(gconj(cstar(N)), rj)
        ok &= is_gint(gsub(v, H))
        for a, b in product((1, 2), (1, 2)):
            q = (Fr(N * a, 3), Fr(N * b, 3))
            v = gmul(gconj(q), rj)
            mi = (Fr(0), Fr(-1))  # -i
            mij = gpow(mi, j) if j >= 0 else gpow(I, -j)  # (-i)^j
            pred = gmul(gmul((Fr((-1) ** k), Fr(0)), mij), (Fr(a, 3), Fr(-b, 3)))
            ok &= is_gint(gsub(v, pred))
            fr = (v[0] - fl(v[0]), v[1] - fl(v[1]))
            ok &= fr[0] in (Fr(1, 3), Fr(2, 3)) and fr[1] in (Fr(1, 3), Fr(2, 3))
print("   conj(c*)rho^j in h + Z[i] and conj(q)rho^j = (1/3)(-1)^k(-i)^j(a - bi) mod Z[i], coordinates in {1/3,2/3}:", ok)
# p_{k-2} = -p_k
ok2 = True
for k in range(2, 10):
    N = 5 ** k
    for q in QN(N):
        for sg in (1, -1):
            pk = gmul(gconj(q), rho_pow(sg * k))
            pk2 = gmul(gconj(q), rho_pow(sg * (k - 2)))
            ok2 &= is_gint(gadd(pk, pk2))
print("   value at rho^(+-(k-2)) = - value at rho^(+-k) mod Z[i] (k = 2..9):", ok2)

print("== Lemma 0: nonzero classes of Lambda_+/-")
for sg in (1, -1):
    pts = [z for z in lattice_pm(sg, Fr(1, 2)) if abs(z[0]) < Fr(1, 2) and abs(z[1]) < Fr(1, 2)]
    print(f"   Lambda_{'+' if sg > 0 else '-'} minimal representatives: {[(str(z[0]), str(z[1])) for z in pts]}")
    print("     each has a coordinate = +-2/5:", all(Fr(2, 5) in (abs(z[0]), abs(z[1])) for z in pts))

print("== Lemma 2: P_1^(1) and the q-hexagon")
P1 = Pk(1, Fr(1))
claimed = sorted({(sx * a, sy * b) for (a, b) in ((Fr(1), Fr(1, 3)), (Fr(1, 3), Fr(1)), (Fr(5, 7), Fr(5, 7)))
                  for sx in (1, -1) for sy in (1, -1)})
print("   P_1^(1) vertices == claimed 12-gon:", P1.V == claimed, "; max |x|^2 =", max(v[0] ** 2 + v[1] ** 2 for v in P1.V), "(10/9 claimed)")
for s in (Fr(1, 6), Fr(11, 56), Fr(1, 5), Fr(491, 2500)):
    Ps = Pk(1, s)
    print(f"   P_1^({s}) = s * P_1^(1):", Ps.V == sorted((s * v[0], s * v[1]) for v in P1.V))
hexq = sorted([(Fr(-5, 4), Fr(0)), (Fr(-5, 7), Fr(-5, 7)), (Fr(0), Fr(-5, 4)), (Fr(1), Fr(-1, 2)), (Fr(1), Fr(1)), (Fr(-1, 2), Fr(1))])
q5 = (Fr(5, 3), Fr(5, 3))
for eps in (Fr(1, 30), Fr(5, 168), Fr(1, 50), Fr(1, 1000), Fr(1, 29), Fr(1, 25)):
    r = Fr(1, 3) - eps
    Y = Yk(q5, 1, r)
    sc = sorted((v[0] / eps, v[1] / eps) for v in Y.V)
    rad = [max(v[0] ** 2 + v[1] ** 2 for v in Yk(q, 1, r).V) / eps ** 2 for q in QN(5)]
    print(f"   eps = {eps}: Y_1((5/3)(1+i)) / eps == claimed hexagon: {sc == hexq}; max|y|^2/eps^2 over the 4 q: {rad}")

print("== Lemma C: exact LP  min s  s.t. x in P_{k-1}^(s), conj(x) rho^{+-k} + nu in Sq_s (nu in Lambda_+- minus Z[i], |nu_e| <= 3/2)")
for k in range(2, 7):
    for sg in (1, -1):
        best = None
        for nu in lattice_pm(sg):
            cons = []
            for j in range(-(k - 1), k):
                for f in funcs(j):
                    cons.append((f[0], f[1], Fr(-1), Fr(0)))
                    cons.append((-f[0], -f[1], Fr(-1), Fr(0)))
            for e, f in enumerate(funcs(sg * k)):
                cons.append((f[0], f[1], Fr(-1), -nu[e]))
                cons.append((-f[0], -f[1], Fr(-1), nu[e]))
            res = lp_min((0, 0, 1), cons)
            if res is not None and (best is None or res[0] < best[0]):
                best = (res[0], nu, res[1])
        print(f"   k = {k}, sign {sg:+d}: min s = {best[0]} (11/56 = {Fr(11, 56)}) at nu = ({best[1][0]}, {best[1][1]}), x = ({best[2][0]}, {best[2][1]})")

print("== Lemma C, only the constraints used in the written proof: w = u rho^-2 in Sq_s, u + nu in Sq_s, nu minimal; u free")
rm2 = rho_pow(-2)
for nu in lattice_pm(1, Fr(1, 2)):
    # variables u1, u2, s ; w = u * rho^-2
    A, B = rm2
    w1 = (A, -B)  # Re(u (A + iB)) = A u1 - B u2
    w2 = (B, A)   # Im = B u1 + A u2
    cons = []
    for f in (w1, w2):
        cons.append((f[0], f[1], Fr(-1), Fr(0)))
        cons.append((-f[0], -f[1], Fr(-1), Fr(0)))
    for e, f in enumerate(((Fr(1), Fr(0)), (Fr(0), Fr(1)))):
        cons.append((f[0], f[1], Fr(-1), -nu[e]))
        cons.append((-f[0], -f[1], Fr(-1), nu[e]))
    res = lp_min((0, 0, 1), cons)
    print(f"   nu = ({nu[0]}, {nu[1]}): min s = {res[0]}")

