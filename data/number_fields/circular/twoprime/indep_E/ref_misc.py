"""Referee E: (1) the order-41 character (12+12i)/41 under <i,rho,sigma>: minimal margin;
(2) sub-windows: G(1,1) with one rotation class (j,l) removed; components at theta=2/7 and exact kappa,
count of non-type components with kappa > 2/7 (these survive just above 2/7)."""
from fractions import Fraction as Fr
from math import floor
from ref_rot import rotations_from_generators, gmul
from ref_vertex import run
from ref_kappa import margins, kappa_primal, kappa_dual

M = 41
i5, i13 = pow(5, -1, M), pow(13, -1, M)
gens = [((3 * i5) % M, (4 * i5) % M), ((5 * i13) % M, (12 * i13) % M), (0, 1)]
for name, G in (("one-prime", [gens[0], gens[2]]), ("two-prime", gens)):
    orb = {(12, 12)}
    st = [(12, 12)]
    while st:
        x, y = st.pop()
        for a, b in G:
            z = ((a * x - b * y) % M, (a * y + b * x) % M)
            if z not in orb:
                orb.add(z); st.append(z)
    print("(1) %s orbit of 12+12i mod 41: size %d, min distance*41 = %d" % (name, len(orb), min(min(x, M - x, y, M - y) for x, y in orb)))

D, rots = rotations_from_generators(1, 1)
classes = [(j, l) for j in (-1, 0, 1) for l in (-1, 0, 1)]
types = [(Fr(65, 2), Fr(65, 2))] + [(Fr(65 * al, 3), Fr(65 * be, 3)) for al in (1, 2) for be in (1, 2)]
theta = Fr(2, 7)
for rem in classes:
    funcs = []
    for (j, l) in classes:
        if (j, l) == rem:
            continue
        z = rots[(0, j, l)]
        funcs.append(z)
        funcs.append(gmul((0, -1), z))
    found, npts, nin = run(2, 7, 65, funcs)
    def iv_of(P):
        out = []
        for (p, q) in funcs:
            v = (p * P[0] + q * P[1]) / 65
            n = floor(v)
            assert theta <= v - n <= 1 - theta
            out.append(n)
        return tuple(out)
    tiv = set(iv_of(P) for P in types)
    surv = 0
    kmax = None
    for iv in found:
        if iv in tiv:
            continue
        ms = margins(funcs, iv, 65)
        kp, _ = kappa_primal(ms)
        kd, _ = kappa_dual(ms)
        assert kp == kd
        if kp > theta:
            surv += 1
            kmax = kp if kmax is None or kp > kmax else kmax
    print("(2) remove class (j,l)=%s: components at 2/7: %d ; non-type components surviving above 2/7: %d ; max kappa among them: %s"
          % (rem, len(found), surv, kmax))
