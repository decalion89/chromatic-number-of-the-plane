"""Meet in the middle, and stop comparing against a constant.

The chain closes when |w_1 + w_2 + w_3|^2 = 1/3.  Normalising w_1 = 1 by the
rotation freedom, that is |1 + a + b|^2 = 1/3, and expanding it into a
condition on (a + abar) + (b + bbar) + (a.bbar + abar.b) invites exactly the
arithmetic slip that cost the last sweep: the target is -8/3, not -4/3, and
the k = 2 control the sweep carried was right either way, so it passed while
the k = 3 search looked for |1 + a + b|^2 = 5/3.

This formulation cannot make that mistake, and is faster by two orders.  Put
s = 1 + a + b.  Then |3s|^2 = 3, and EVERY element of modulus sqrt3 is
u(1 + zeta_6) for a unit step u, because (1 + zeta_6)(1 + zeta_6bar) = 3 and
the quotient has modulus one.  So the reachable s are

    s = u(1 + zeta_6)/3,     u a unit step,

one field multiplication each.  Then a + b = s - 1 is known, and the search is
a hash lookup: for every a, ask whether s - 1 - a is itself a unit step.  No
target constant appears anywhere.
"""
import sys, itertools, time, pickle
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

K = CycloField(33)
D = K.degree
ONE = K.rational(1)


def inv(a):
    rows = [list(K.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        s = Fraction(1) / M[c][c]
        M[c] = [v * s for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


Z11, Z3 = K.zeta(3), K.zeta(11)
Z6 = K.neg(K.mul(Z3, Z3))
g, p = K.zero(), ONE
for k in range(1, 11):
    p = K.mul(p, Z11) if k > 1 else Z11
    g = K.add(g, p if k in {1, 3, 4, 5, 9} else K.neg(p))
assert K.mul(g, g) == K.rational(-11)
RHO = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))
assert K.norm2(RHO) == ONE
ONE_PLUS = K.add(ONE, Z6)
assert K.norm2(ONE_PLUS) == K.rational(3)

t0 = time.time()
extra = set()
for coeffs in itertools.product(range(-1, 2), repeat=6):
    if not any(coeffs):
        continue
    a, q = K.zero(), ONE
    for c in coeffs:
        if c:
            a = K.add(a, tuple(Fraction(c) * x for x in q))
        q = K.mul(q, Z11)
    try:
        u = K.mul(a, inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) == ONE:
        extra.add(u)
zp, z = [], ONE
for _ in range(6):
    zp.append(z)
    z = K.mul(z, Z6)
steps = set()
for c in zp:
    for e in {ONE} | extra:
        base = K.mul(c, e)
        steps.add(base)
        steps.add(K.mul(base, RHO))
        steps.add(K.mul(base, K.conj(RHO)))
        steps.add(K.mul(base, K.mul(RHO, RHO)))
        steps.add(K.mul(base, K.mul(K.conj(RHO), K.conj(RHO))))
steps = sorted(steps)
sset = set(steps)
print(f"{len(steps)} unit steps  [{time.time()-t0:.0f}s]", flush=True)

# control: k = 2 closes iff |1 + a|^2 = 1/3, i.e. a = -(5 -+ sqrt-11)/6
ctrl = [a for a in steps if K.norm2(K.add(ONE, a)) == K.rational(Fraction(1, 3))]
print(f"k = 2 control: {len(ctrl)} Moser rotations", flush=True)

THIRD = Fraction(1, 3)
found = []
for i, u in enumerate(steps):
    s = tuple(THIRD * x for x in K.mul(u, ONE_PLUS))
    z = K.sub(s, ONE)
    for a in steps:
        b = K.sub(z, a)
        if b in sset:
            found.append((a, b))
    if i % 400 == 0:
        print(f"  ... {i}/{len(steps)}, {len(found)} chains "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(found)} ordered k = 3 chains  [{time.time()-t0:.0f}s]", flush=True)
with open(f"{sys.path[0]}/chains3.pkl", "wb") as fh:
    pickle.dump([(list(a), list(b)) for a, b in found], fh)
