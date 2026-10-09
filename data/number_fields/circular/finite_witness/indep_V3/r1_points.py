"""Referee check 1/2: exact arithmetic on the points of the two nine-point witnesses.

Elements of Q(sqrt d) are represented as pairs (x, y) of Fractions meaning x + y*sqrt(d).  Nothing is imported from
the repository; the witness files are read as JSON only.
"""
import gzip, json, sys, itertools
from fractions import Fraction as Fr

BASE = sys.argv[1]  # folder holding the witness files


class QF:
    """x + y sqrt(d), exact."""
    def __init__(self, x, y, d):
        self.x, self.y, self.d = Fr(x), Fr(y), d
    def __sub__(s, o): return QF(s.x - o.x, s.y - o.y, s.d)
    def __add__(s, o): return QF(s.x + o.x, s.y + o.y, s.d)
    def __mul__(s, o): return QF(s.x * o.x + s.d * s.y * o.y, s.x * o.y + s.y * o.x, s.d)
    def __eq__(s, o): return s.x == o.x and s.y == o.y
    def __hash__(s): return hash((s.x, s.y))
    def __repr__(s): return f'({s.x} + {s.y} r)'


def load(name):
    with gzip.open(f'{BASE}/{name}', 'rt') as f:
        return json.load(f)


def points_of(W):
    d, D = W['d'], W['denominator']
    return [(QF(Fr(a, D), Fr(b, D), d), QF(Fr(c, D), Fr(e, D), d)) for a, b, c, e in W['points']]


def dist2(P, Q):
    dx, dy = P[0] - Q[0], P[1] - Q[1]
    return dx * dx + dy * dy


def is_squarefree(d):
    return all(d % (k * k) for k in range(2, int(d ** 0.5) + 2))


ONE = lambda d: QF(1, 0, d)

for name in ['witness_q7.json.gz', 'witness_q31.json.gz']:
    W = load(name)
    d = W['d']
    print(f'== {name}: d = {d}, D = {W["denominator"]}, p/q = {W["p"]}/{W["q"]}, keys = {sorted(W)}')
    print('   d squarefree:', is_squarefree(d), ' d not a square:', int(d ** 0.5) ** 2 != d)
    P = points_of(W)
    n = len(P)
    print('   n =', n, ' distinct:', len(set(P)) == n)
    units = sorted((i, j) for i, j in itertools.combinations(range(n), 2) if dist2(P[i], P[j]) == ONE(d))
    E = sorted(tuple(sorted(e)) for e in W['edges'])
    print('   listed edges:', len(E), ' unit pairs:', len(units), ' equal (induced, all unit):', units == E)
    print('   no repeated edge:', len(set(E)) == len(E))
    # all pairwise squared distances, to see none is accidentally 1
    for i, j in itertools.combinations(range(n), 2):
        dd = dist2(P[i], P[j])
        if (i, j) not in units and dd == ONE(d):
            print('   MISSING unit pair', i, j)
    # check no unit triangle within the set (sanity)
    tri = [t for t in itertools.combinations(range(n), 3)
           if all(tuple(sorted(e)) in set(units) for e in itertools.combinations(t, 2))]
    print('   triangles among unit pairs:', tri)
    for k, p in enumerate(P):
        print(f'   v{k}: x = {p[0].x} + {p[0].y} r, y = {p[1].x} + {p[1].y} r')

# Paper coordinates of P0..P7, S, typed by hand from the tex source (x = a + b r, y = c + e r, r = sqrt 7)
paper = [
    ((-1, 0), (0, 0)),                                   # P0 = (-1, 0)
    ((Fr(-3, 8), Fr(-1, 8)), (Fr(-5, 8), Fr(-1, 8))),     # P1 = (-(3+r)/8, -(5+r)/8)
    ((Fr(1, 4), 0), (0, Fr(-1, 4))),                     # P2 = (1/4, -r/4)
    ((0, Fr(1, 4)), (Fr(1, 4), 0)),                      # P3 = (r/4, 1/4)
    ((Fr(-1, 4), 0), (0, Fr(1, 4))),                     # P4 = (-1/4, r/4)
    ((Fr(3, 8), Fr(-1, 8)), (Fr(-5, 8), Fr(1, 8))),      # P5 = ((3-r)/8, (r-5)/8)
    ((1, 0), (0, 0)),                                    # P6 = (1, 0)
    ((0, 0), (0, 0)),                                    # P7 = (0, 0)
    ((0, 0), (1, 0)),                                    # S  = (0, 1)
]
W = load('witness_q7.json.gz'); P = points_of(W)
ok = all(P[k][0] == QF(*paper[k][0], 7) and P[k][1] == QF(*paper[k][1], 7) for k in range(9))
print('paper coordinates P0..P7, S equal to the file (in order 0..8):', ok)

# README coordinates of the Q(sqrt31) witness, in file order 0..8
readme31 = [
    ((Fr(-8, 5), 0), (Fr(-4, 5), 0)),
    ((-1, 0), (0, 0)),
    ((Fr(-60, 80), Fr(-3, 80)), (Fr(45, 80), Fr(-4, 80))),
    ((Fr(-3, 5), 0), (Fr(-4, 5), 0)),
    ((0, 0), (0, 0)),
    ((0, Fr(1, 16)), (Fr(-15, 16), 0)),
    ((Fr(-30, 40), Fr(1, 40)), (Fr(-15, 40), Fr(-2, 40))),
    ((Fr(-108, 80), Fr(-3, 80)), (Fr(-19, 80), Fr(-4, 80))),
    ((-1, Fr(1, 16)), (Fr(-15, 16), 0)),
]
W = load('witness_q31.json.gz'); P = points_of(W)
ok = all(P[k][0] == QF(*readme31[k][0], 31) and P[k][1] == QF(*readme31[k][1], 31) for k in range(9))
print('README coordinates of the Q(sqrt31) witness equal to the file:', ok)
