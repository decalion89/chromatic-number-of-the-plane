"""Referee E: spot checks of Remark (3) of the note with the vertex method (ref_vertex.py) and exact kappa
(ref_kappa.py primal/dual LP).  Usage: python3 ref_remark3.py WINDOW NUM DEN [kappa]"""
import sys
from fractions import Fraction as Fr
from math import floor
from ref_vertex import functionals_for, run
from ref_kappa import margins, kappa_primal, kappa_dual

window, a, b = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
dokappa = len(sys.argv) > 4
D, funcs = functionals_for(window)
found, npts, nin = run(a, b, D, funcs)
theta = Fr(a, b)
# type points at modulus D
def iv_of(P):
    out = []
    for (p, q) in funcs:
        v = (p * P[0] + q * P[1]) / D
        n = floor(v)
        if not (theta <= v - n <= 1 - theta):
            return None
        out.append(n)
    return tuple(out)
types = [(Fr(D, 2), Fr(D, 2))] + [(Fr(D * al, 3) % D, Fr(D * be, 3) % D) for al in (1, 2) for be in (1, 2)]
tiv = set(iv_of(P) for P in types)
assert None not in tiv
ntype = sum(1 for iv in found if iv in tiv)
print("window %s (D=%d, %d functionals), theta=%s: components %d, of which contain a type point: %d, extra: %d"
      % (window, D, len(funcs), theta, len(found), ntype, len(found) - ntype))
if dokappa:
    ks = {}
    for iv in found:
        if iv in tiv:
            continue
        ms = margins(funcs, iv, D)
        kp, _ = kappa_primal(ms)
        kd, _ = kappa_dual(ms)
        assert kp == kd
        ks[kp] = ks.get(kp, 0) + 1
    print("  exact kappa of the extra components (value: count):", {str(k): v for k, v in sorted(ks.items())})
    print("  max kappa among extras:", max(ks) if ks else None)
