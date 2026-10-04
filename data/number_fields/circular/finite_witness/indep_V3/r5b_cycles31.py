"""Count 5-cycles 0 - a - b - c - d - 0 over Q(sqrt31) with all points of denominator 80, a, d unit vectors,
b, c not in {0} u U (two extra points), up to the symmetries of {0} u U (x<->y, sign changes, sqrt31 -> -sqrt31,
and the reversal of the cycle)."""
import json, itertools
D, d = 80, 31
U = [tuple(u) for u in json.load(open('seed31_cands.json'))['units']]
Us = set(U)
add = lambda p, q: tuple(x + y for x, y in zip(p, q))
zero = (0, 0, 0, 0)
# b = a + u (b not in U, b != 0)
nb = {}
for a in U:
    for u in U:
        b = add(a, u)
        if b != zero and b not in Us:
            nb.setdefault(b, set()).add(a)
cyc = set()
for b, As in nb.items():
    for u in U:
        c = add(b, u)
        if c in nb and c != b:
            for a in As:
                for dd in nb[c]:
                    if a != dd:
                        cyc.add((a, b, c, dd))
def syms(p):
    a, b, c, e = p
    out = []
    for (x, y) in [((a, b), (c, e)), ((c, e), (a, b))]:
        for sx in (1, -1):
            for sy in (1, -1):
                for g in (1, -1):
                    out.append((sx * x[0], sx * g * x[1], sy * y[0], sy * g * y[1]))
    return out
def canon(t):
    best = None
    for k in range(8):
        im = [syms(p)[k] for p in t]
        for seq in (tuple(im), tuple(reversed(im))):
            if best is None or seq < best:
                best = seq
    return best
classes = {}
for t in cyc:
    classes.setdefault(canon(t), t)
print('directed 5-cycles through 0 with two extra points:', len(cyc), '; classes up to symmetry and reversal:', len(classes))
json.dump([list(map(list, t)) for t in classes.values()], open('cycles31_classes.json', 'w'))
