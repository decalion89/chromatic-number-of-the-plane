"""E(t, L, B): 'there is an independent dominating set S of G_13 with |S| >= t (and <= tmax) whose indicator is
lex-leader under Aut(G_13) on the first L positions of lex_order(), and S is none of the blocked sets'.

Blocking: for every listed set R (any member of its orbit), every image of R whose L-prefix is the largest in its
orbit (exactly the images that can satisfy the lex-leader clauses) gets the clause OR_{v in image} -s_v.  Since S
is independent and dominating (= maximal independent), S containing an image forces S = image.  So E is
unsatisfiable iff every maximal independent set with t..tmax points lies in the orbit of a listed set.
Variables: s_v = 1 + v; the rest auxiliary.

usage: enum_cert.py OUT.cnf t [--tmax T] [--L 25] [--sets FILE ...]
"""
import argparse
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from g13 import NV, EDGES, ADJ, POINTS, INDEX, add, circle, automorphisms
from g13cnf import CNF, totalizer_atleast, totalizer_atmost, lex_order, lex_leader_sets

AUT = np.array(automorphisms(), dtype=np.int32)          # AUT[g][v] = image of v


def leader_images(S, order, L):
    """images g(S) whose indicator, read on order[:L], is maximal over the orbit"""
    ind = np.zeros((len(AUT), NV), dtype=np.int8)
    rows = np.repeat(np.arange(len(AUT)), len(S))
    ind[rows, AUT[:, S].ravel()] = 1
    pref = ind[:, order[:L]]
    # lexicographic max of rows (1 > 0): sort descending
    keys = [pref[:, i] for i in range(L - 1, -1, -1)]
    idx = np.lexsort(keys)
    best = pref[idx[-1]]
    hit = np.where((pref == best).all(axis=1))[0]
    out = {frozenset(np.nonzero(ind[g])[0].tolist()) for g in hit}
    return [sorted(x) for x in out]


ROSETTE_CIRCLES = (6, 7, 9, 11)       # the circles N = c other than N = 1 with no unit distance inside


def build(t, tmax, L, sets, norosette=False, rosette=None):
    cnf = CNF()
    cnf.top = NV
    s = lambda v: 1 + v
    for a, b in EDGES:
        cnf.add([-s(a), -s(b)])
    for v in range(NV):
        cnf.add([s(v)] + [s(w) for w in ADJ[v]])
    totalizer_atleast(cnf, [s(v) for v in range(NV)], t)
    if tmax:
        totalizer_atmost(cnf, [s(v) for v in range(NV)], tmax)
    order = lex_order()
    perms = [list(p) for p in AUT if any(p[v] != v for v in range(NV))]
    if L:
        lex_leader_sets(cnf, {v: s(v) for v in range(NV)}, perms, order, L)
    if norosette:
        # Part B: no point of S together with its whole circle N = c, c in ROSETTE_CIRCLES
        for c in ROSETTE_CIRCLES:
            cc = circle(c)
            for p in POINTS:
                cnf.add([-s(INDEX[p])] + [-s(INDEX[add(p, q)]) for q in cc])
    if rosette:
        # Part A_c: 0 and the whole circle N = c lie in S (lex-leader clauses must then be dropped: L = 0)
        for q in [(0, 0)] + circle(rosette):
            cnf.add([s(INDEX[q])])
    nblock = 0
    for S in sets:
        for img in leader_images(S, order, L):
            cnf.add([-s(v) for v in img])
            nblock += 1
    return cnf, nblock


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("t", type=int)
    ap.add_argument("--tmax", type=int, default=0)
    ap.add_argument("--L", type=int, default=25)
    ap.add_argument("--sets", nargs="*", default=[])
    ap.add_argument("--norosette", action="store_true", help="Part B: forbid a point with its whole circle 6/7/9/11")
    ap.add_argument("--rosette", type=int, default=0, help="Part A_c: 0 and circle N = c in S (use --L 0)")
    a = ap.parse_args()
    sets = []
    for fn in a.sets:
        sets += [list(map(int, l.split())) for l in open(fn) if l.strip()]
    sets = [S for S in sets if len(S) >= a.t and (not a.tmax or len(S) <= a.tmax)]
    cnf, nb = build(a.t, a.tmax, a.L, sets, a.norosette, a.rosette)
    with open(a.out, "w") as fh:
        fh.write(cnf.text())
    print(f"{a.out}: t={a.t} tmax={a.tmax} L={a.L}: {len(sets)} listed sets, {nb} blocking clauses; "
          f"{cnf.top} vars, {len(cnf.clauses)} clauses")


if __name__ == "__main__":
    main()
