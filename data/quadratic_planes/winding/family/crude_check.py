"""Exact checks of the elementary facts behind the crude form of Proposition 1 (the form proved in
lean/FourColours.lean; notes/four_colours_11_mod_12.md, section 7):

* (3 + 4i)^m = A + Bi has A odd and B even; (3 + 4i)^m = 3 - i (mod 5) and (3 - 4i)^m = 3 + i (mod 5) for m >= 1;
  (3 + 4i)^m = i^m (mod 3);
* with t = p + qi and conj(t)(3 - i) = a + bi: a = 2b (mod 5), and a = b = 0 (mod 5) iff q = 3p (mod 5); with
  conj(t)(3 + i) = a + bi: a = -2b (mod 5), and a = b = 0 (mod 5) iff q = 2p (mod 5);
* type q: if gamma_1/3 + a/5 and gamma_2/3 + b/5 both lie in [1/3, 2/3] + Z, with gamma_j in {1, 2} and
  a = +-2b (mod 5), then a = b = 0 (mod 5);
* type c: if a = +-2b (mod 5) and b != 0 (mod 5), then (s1 - a/5)^2 + (s2 - b/5)^2 >= 1/18 for all
  s1, s2 in [-1/6, 1/6] (equality at a = 2, b = 1).

Run: python3 crude_check.py   (prints ALL PLAN CHECKS PASSED)
"""
from fractions import Fraction as F
def mul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
z = (1, 0)
for m in range(1, 60):
    z = mul(z, (3, 4))
    assert (z[0] % 5, z[1] % 5) == (3, 4), m                       # (3+4i)^m = 3 - i (mod 5)
    im = [(1,0),(0,1),(-1,0),(0,-1)][m % 4]
    assert ((z[0]-im[0]) % 3, (z[1]-im[1]) % 3) == (0, 0), m        # (3+4i)^m = i^m (mod 3)
    assert z[0] % 2 == 1 and z[1] % 2 == 0, m                      # A odd, B even
    zb = (z[0], -z[1])
    assert (zb[0] % 5, zb[1] % 5) == (3, 1), m                      # (3-4i)^m = 3 + i (mod 5)
# t = p + q i, conj t = p - q i; (3 - i)(p - q i) = a + b i: a = 2b (mod 5); zero iff q = 3p (mod 5)
for p in range(5):
    for q in range(5):
        a, b = mul((3, -1), (p, -q)); assert (a - 2*b) % 5 == 0
        assert ((a % 5, b % 5) == (0, 0)) == (q % 5 == (3*p) % 5)
        a2, b2 = mul((3, 1), (p, -q)); assert (a2 + 2*b2) % 5 == 0
        assert ((a2 % 5, b2 % 5) == (0, 0)) == (q % 5 == (2*p) % 5)
# type q residue lemma: gamma/3 + a/5 in [1/3, 2/3] + Z
def ok(g, a): return 5 <= (5*g + 3*a) % 15 <= 10
for g1 in (1, 2):
    for g2 in (1, 2):
        for b in range(5):
            for sgn in (2, -2):
                a = (sgn*b) % 5
                if ok(g1, a) and ok(g2, b): assert a == 0 and b == 0, (g1, g2, a, b)
# distance lemma: a = +-2b (mod 5), b != 0 (mod 5), s in [-1/6,1/6]^2  =>  (s1 - a/5)^2 + (s2 - b/5)^2 >= 1/18
def d1(t):  # distance from t to [-1/6, 1/6]
    return max(F(0), abs(t) - F(1, 6))
mn = None
for a in range(-30, 31):
    for b in range(-30, 31):
        if b % 5 == 0 or ((a - 2*b) % 5 and (a + 2*b) % 5): continue
        v = d1(F(a, 5))**2 + d1(F(b, 5))**2
        assert v >= F(1, 18), (a, b, v)
        mn = v if mn is None or v < mn else mn
print("min distance^2 =", mn, "(1/18 =", F(1, 18), ")")
# base case: x in [-1/6,1/6]^2 with |x|^2 = 1/18 only at corners
print("ALL PLAN CHECKS PASSED")
