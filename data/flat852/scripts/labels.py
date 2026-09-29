"""labels.py -- analyse a 4-colouring: edge label = the partition {{c(x),c(y)}, complement} of the 4 colours
(3 possible labels; equal labels <=> equal-or-disjoint colour pairs; invariant under colour permutations).
A colouring is linear on a part iff all edges of each direction (up to sign) carry the same label there."""
import sys, json
import numpy as np
from collections import Counter, defaultdict
from flat import *
from cosets import cls

LAB = {}
for a in range(4):
    for b in range(4):
        if a != b:
            p = frozenset((a, b))
            LAB[(a, b)] = min(p, frozenset(range(4)) - p, key=lambda s: sorted(s))
LNAME = {frozenset((0, 1)): "A", frozenset((0, 2)): "B", frozenset((0, 3)): "C"}


def analyse(P, col, U, part_rows=None):
    E, J = build_edges(P, U)
    C = cls(P)
    by = defaultdict(Counter)       # (part, direction up to sign) -> label counts
    neg = {}
    for j, u in enumerate(U):
        neg[j] = U.index(tuple(-x for x in u))
    for (a, b), j in zip(E, J):
        ca, cb = C[a], C[b]
        jj = min(j, neg[j])
        part = ("row", ca[1]) if ca[1] == cb[1] else ("col", ca[0])
        if ca == cb:
            part = ("cell", ca)
        by[(part, jj)][LNAME[LAB[(col[a], col[b])]]] += 1
    return by


if __name__ == "__main__":
    P = np.load(sys.argv[1]); st = json.load(open(sys.argv[2]))
    col = st["col"]
    D, CD, U = directions()
    P = P[:len(col)]
    by = analyse(P, col, U)
    parts = sorted(set(k[0] for k in by), key=str)
    for part in parts:
        items = [(j, c) for (p, j), c in by.items() if p == part]
        tot = sum(sum(c.values()) for j, c in items)
        mixed = [(j, dict(c)) for j, c in items if len(c) > 1]
        nmin = sum(sum(c.values()) - max(c.values()) for j, c in items)
        print(f"{part}: {len(items)} directions, {tot} edges; directions with mixed labels: {len(mixed)}; "
              f"minority-label edges: {nmin}")
