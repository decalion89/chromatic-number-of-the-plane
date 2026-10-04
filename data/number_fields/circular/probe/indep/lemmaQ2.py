"""Lemma Q exactly as written: hypothesis |y| < 1/6 is used only through |u_e| <= |u| = |y| < 1/6 (u = conj(y) rho^k),
so relax it to the box |u_e| <= 1/6 on u (implied by |y| < 1/6).  Cells of q at rho^(k-2) (two-sided), condition at
rho^k with nu' over all of Lambda_+ minus Z[i] in a box.  Also print the eps = 0 witness when the box is put on y."""
from fractions import Fraction as Fr
from hgeom import funcs, ev, fl, gmul, gconj, rho_pow
from lp3x import lp_min
from sn import QN
from lemmas_a import lattice_pm


def lp(q, k, sg, boxon):
    best = None
    pk = gmul(gconj(q), rho_pow(sg * k))
    for nu in lattice_pm(sg, Fr(1)):
        cons = [(Fr(0), Fr(0), Fr(-1), Fr(0))]
        boxf = funcs(sg * k) if boxon == "u" else [(Fr(1), Fr(0)), (Fr(0), Fr(1))]
        for f in boxf:
            cons.append((f[0], f[1], Fr(0), Fr(1, 6)))
            cons.append((-f[0], -f[1], Fr(0), Fr(1, 6)))
        for f in funcs(sg * (k - 2)):
            val = ev(f, q); n = fl(val)
            cons.append((f[0], f[1], Fr(-1), n + Fr(2, 3) - val))
            cons.append((-f[0], -f[1], Fr(-1), val - n - Fr(1, 3)))
        for e, f in enumerate(funcs(sg * k)):
            base = (pk[e] - fl(pk[e])) + nu[e]
            cons.append((f[0], f[1], Fr(-1), Fr(2, 3) - base))
            cons.append((-f[0], -f[1], Fr(-1), base - Fr(1, 3)))
        res = lp_min((0, 0, 1), cons)
        if res is not None and (best is None or res[0] < best[0]):
            best = (res[0], nu, res[1])
    return best


for k in (2, 3, 4, 5):
    N = 5 ** k
    vals = set()
    for q in QN(N):
        for sg in (1, -1):
            vals.add(lp(q, k, sg, "u")[0])
    print(f"k = {k}: box on u: min eps over 4 q, both signs = {sorted(vals)}")
q = QN(25)[0]
b = lp(q, 2, 1, "y")
y = b[2]
print("box on y instead (k=2, q=25/3(1+i), sign +): min eps =", b[0], "at y =", (str(y[0]), str(y[1])), "|y|^2 =", y[0]**2 + y[1]**2, "(1/36 =", Fr(1,36), ") nu' =", (str(b[1][0]), str(b[1][1])))
