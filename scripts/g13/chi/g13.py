"""G_13 = Cay(F_169, mu_14): vertices z = (x, y) = x + y w, w^2 = 2, index v = 13 x + y; z ~ z' iff N(z - z') = 1,
N(x, y) = x^2 - 2 y^2 mod 13.  Shared helpers for the 5-colouring experiments."""
Q, N0 = 13, 2
K = 5


def norm(z):
    return (z[0] * z[0] - N0 * z[1] * z[1]) % Q


def add(a, b):
    return ((a[0] + b[0]) % Q, (a[1] + b[1]) % Q)


def sub(a, b):
    return ((a[0] - b[0]) % Q, (a[1] - b[1]) % Q)


def neg(a):
    return ((-a[0]) % Q, (-a[1]) % Q)


def mul(a, b):
    return ((a[0] * b[0] + N0 * a[1] * b[1]) % Q, (a[0] * b[1] + a[1] * b[0]) % Q)


def conj(a):
    return (a[0], (-a[1]) % Q)


POINTS = [(x, y) for x in range(Q) for y in range(Q)]
INDEX = {z: Q * z[0] + z[1] for z in POINTS}
NV = len(POINTS)
UNITS = [z for z in POINTS if norm(z) == 1]
EDGES = sorted({tuple(sorted((INDEX[z], INDEX[add(z, u)]))) for z in POINTS for u in UNITS})
ADJ = [sorted(INDEX[add(z, u)] for u in UNITS) for z in POINTS]


def order(z):
    k, w = 1, z
    while w != (1, 0):
        w, k = mul(w, z), k + 1
    return k


GEN = next(u for u in UNITS if order(u) == Q + 1)          # generator of mu_14
POW = [(1, 0)]
for _ in range(Q):
    POW.append(mul(POW[-1], GEN))                          # POW[k] = GEN^k, k = 0..13


def circle(c):
    return [z for z in POINTS if norm(z) == c]


def stabiliser():
    """the 28 maps fixing 0: z -> l z (k = 0..13, l = GEN^k) then z -> l conj(z); as vertex permutations"""
    out = []
    for c in (0, 1):
        for l in POW:
            out.append([INDEX[mul(l, conj(z) if c else z)] for z in POINTS])
    return out


def translation(a):
    return [INDEX[add(z, a)] for z in POINTS]


def automorphisms():
    """all 169 * 28 maps z -> s(z) + a"""
    st = stabiliser()
    return [[INDEX[add(POINTS[p[v]], a)] for v in range(NV)] for a in POINTS for p in st]


def check_colouring(col, k=K):
    assert len(col) == NV and all(0 <= c < k for c in col)
    bad = [(a, b) for a, b in EDGES if col[a] == col[b]]
    return bad
