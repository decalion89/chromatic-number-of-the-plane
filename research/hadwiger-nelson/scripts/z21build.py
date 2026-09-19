"""Q(zeta_21): triangles from zeta_6 and blocking inherited from zeta_7.

Blocking passes upward through a field extension for free.  If phi were a
homomorphism on the big module, nonzero on every edge vector, its restriction
to the Q(zeta_7) submodule would be one there, nonzero on the 87 directions --
and none exists.  So no sampling is needed to know Q(zeta_21) blocks: it
contains Q(zeta_7), and blocking is inherited.

What zeta_21 adds is zeta_6, hence the Eisenstein lattice, hence triangles --
the thing Q(zeta_7) provably cannot have.  The step set is therefore the union
of the two: sixth roots of unity to fold and to triangulate, Q(zeta_7)'s
modulus-one steps to block.

Both halves are checked here rather than assumed, and the graph is handed to a
solver.
"""
import sys, itertools, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

F = CycloField(21)
D = F.degree
one = F.rational(1)
print(f"Q(zeta_21): degree {D}", flush=True)

z3 = F.zeta(7)                      # zeta_21^7 = zeta_3
z6 = F.neg(F.mul(z3, z3))           # -zeta_3^2 is a primitive sixth root
assert F.add(F.sub(F.mul(z6, z6), z6), one) == F.zero(), "z6^2 - z6 + 1 = 0"
assert F.norm2(z6) == one
print("  zeta_6 present: Eisenstein triangles available", flush=True)

# zeta_7 sits inside as zeta_21^3; its modulus-one steps embed directly.
z7 = F.zeta(3)
assert F.norm2(z7) == one


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def denom(u):
    d = 1
    for x in u:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


# The Eisenstein half: powers of zeta_6.
eis, z = [], one
for _ in range(6):
    z = F.mul(z, z6)
    eis.append(z)

# The blocking half: alpha / conj(alpha) with alpha in Z[zeta_7], which is
# where denominator-29 modulus-one elements come from.
blk = set()
for coeffs in itertools.product(range(-3, 4), repeat=4):
    if not any(coeffs):
        continue
    a = F.zero()
    p = one
    for c in coeffs:
        a = F.add(a, tuple(Fraction(c) * x for x in p))
        p = F.mul(p, z7)
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) != one:
        continue
    d = denom(u)
    if d % 5 and d <= 29:
        blk.add(u)
        blk.add(F.neg(u))
steps = sorted(set(eis) | blk)
print(f"  {len(eis)} Eisenstein steps + {len(blk)} blocking steps "
      f"= {len(steps)} total", flush=True)

pts = {F.zero()}
for rnd in range(2):
    pts = {F.add(p, s) for p in pts for s in steps} | pts
    print(f"    round {rnd+1}: {len(pts)} points", flush=True)

P = sorted(pts)
idx = {p: i for i, p in enumerate(P)}
edges = sorted({(min(i, j), max(i, j))
                for i, p in enumerate(P) for s in steps
                for j in [idx.get(F.add(p, s))] if j is not None and j != i})
print(f"  graph: {len(P)} vertices, {len(edges)} edges", flush=True)

adj = [set() for _ in P]
for a, b in edges:
    adj[a].add(b)
    adj[b].add(a)
tri = sum(len(adj[a] & adj[b]) for a, b in edges) // 3
print(f"  triangles: {tri}", flush=True)

ev = {tuple(F.sub(P[j], P[i])) for i, j in edges}
den = 1
for v in ev:
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
iv = set()
for v in ev:
    w = tuple(int(x * den) for x in v)
    g = 0
    for t in w:
        g = gcd(g, abs(t))
    iv.add(tuple(t // g for t in w) if g > 1 else w)
phi, why = has_homomorphism(sorted(iv), 5)
print(f"  {len(iv)} edge vectors; coset 5-colouring: "
      + ("EXISTS" if phi else "NONE -- blocked"), flush=True)

for k in (3, 4, 5, 6):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  {k}-colourable: {r}" + ("" if r else f"   *** chi > {k} ***"),
          flush=True)
    if r:
        break
