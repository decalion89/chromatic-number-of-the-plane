"""Referee E: the 8 extra components are the 7-adic characters of Prop. C1.
For c = 65(a+bi)/7 and gamma = (p+qi)/65 in G(1,1): Re(conj(c) gamma) = Re((a-bi)(p+qi))/7 = (ap+bq)/7 mod 1.
A 7-adic character: z -> lam(z mod 7)/7 with lam F_7-linear on F_49 = F_7(i), lam(mu_8) in {2,3,4,5}.
Check: (1) each of the 8 points gives such a lam (values on gamma in {2,..,5}/7, so min margin 2/7);
(2) lam(x+yi) = 2x+3y corresponds to 65(1+5i)/7; (3) the F_7-linear maps with lam(mu_8) in {2..5} are exactly 8."""
from fractions import Fraction as Fr
from ref_rot import rotations_from_generators
D, rots = rotations_from_generators(1, 1)
inv65 = pow(65, -1, 7)
mu8 = [(x, y) for x in range(7) for y in range(7) if (x * x + y * y) % 7 == 1]
assert len(mu8) == 8
good = [(u, v) for u in range(7) for v in range(7) if all((u * x + v * y) % 7 in (2, 3, 4, 5) for (x, y) in mu8)]
print("F_7-linear lam(x+yi)=ux+vy with lam(mu_8) in {2..5}:", good)
pts = [(1, 2), (1, 5), (2, 1), (2, 6), (5, 1), (5, 6), (6, 2), (6, 5)]
for (a, b) in pts:
    vals = set()
    for (p, q) in rots.values():
        vals.add((a * p + b * q) % 7)
    # matching lam: Re(conj(c) gamma) = (a p + b q)/7 and gamma mod 7 = (p + q i) * 65^{-1} -> lam(x+yi) = 65(a x + b y)
    lam = ((65 * a) % 7, (65 * b) % 7)
    print("c=65(%d+%di)/7: values 7*Re(conj(c)gamma) mod 7 over G(1,1): %s ; corresponds to lam(x+yi)=%dx+%dy ; in list: %s"
          % (a, b, sorted(vals), lam[0], lam[1], lam in good))
