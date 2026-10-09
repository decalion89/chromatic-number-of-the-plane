"""Referee: unit vectors with a given denominator (own brute force), the q7 seed 5-cycle, and candidate seeds for
Q(sqrt31) with denominator 80 built like q7_seed.json (0, all unit vectors, and two points closing a 5-cycle
0, (-1,0), (-2,0), Q, R with R a unit vector)."""
import json, math, sys, itertools


def units(d, D):
    out = []
    B = math.isqrt(D * D // d)
    for b in range(-B, B + 1):
        for e in range(-B, B + 1):
            r = D * D - d * (b * b + e * e)
            if r < 0:
                continue
            for a in range(-D, D + 1):
                c2 = r - a * a
                if c2 < 0:
                    continue
                c = math.isqrt(c2)
                if c * c != c2:
                    continue
                for cc in {c, -c}:
                    if a * b + cc * e == 0:
                        out.append((a, b, cc, e))
    return sorted(set(out))


def isunit(p, q, d, D):
    a, b, c, e = (x - y for x, y in zip(p, q))
    return a * a + d * b * b + c * c + d * e * e == D * D and a * b + c * e == 0


for d, D in [(7, 160), (31, 80)]:
    U = units(d, D)
    print(f'Q(sqrt{d}), D = {D}: {len(U)} unit vectors; irrational ones: {sum(1 for u in U if u[1] or u[3])}')

# q7 seed 5-cycle
d, D = 7, 160
cyc = [(0, 0, 0, 0), (-160, 0, 0, 0), (-320, 0, 0, 0), (-220, -20, -100, -20), (-120, 0, 0, -40)]
print('q7 seed 5-cycle closes:', all(isunit(cyc[i], cyc[(i + 1) % 5], d, D) for i in range(5)),
      '; (-3/4, -sqrt7/4) = (-120, 0, 0, -40)/160')
U7 = set(units(7, 160))
seed = json.load(open(sys.argv[1]))
S = {tuple(p) for p in seed['points']}
print('q7_seed.json = {0} + all unit vectors + the two extra points:',
      S == U7 | {(0, 0, 0, 0), (-320, 0, 0, 0), (-220, -20, -100, -20)})

# candidates for Q(sqrt31)
d, D = 31, 80
U = units(d, D)
m2 = (-2 * D, 0, 0, 0)
cands = []
for R in U:
    for u in U:
        Q = tuple(x + y for x, y in zip(m2, u))
        if Q in set(U) or Q == (0, 0, 0, 0) or Q == (-D, 0, 0, 0):
            continue
        if isunit(Q, R, d, D) and R != (-D, 0, 0, 0):
            # 5-cycle 0, (-1,0), (-2,0), Q, R; induced? (no chords: automatically, as no unit triangles)
            cands.append((Q, R))
print(f'5-cycles 0, (-1,0), (-2,0), Q, R over Q(sqrt31) with denominator 80: {len(cands)}')
for Q, R in cands[:12]:
    print('   Q =', Q, ' R =', R)
json.dump({'units': U, 'cands': cands}, open('seed31_cands.json', 'w'))
