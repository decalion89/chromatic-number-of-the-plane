"""linstart.py -- linear colourings as warm starts.
pi(x) = 7x mod 2 (12 bits) is an isomorphism M/2M -> O/2O. For a pair (a,b) of 12-bit rows, the colouring
c(x) = 2*(a.pi(x)) + (b.pi(x)) is proper on every edge whose residue r has (a.r, b.r) != 0.
D-proper pairs (42 row spaces): proper on all D-edges (all of L's structure); conj-proper: on all conj(D)-edges.
Mixed start: the D-linear colouring on the L-row of the triangle, the conj-linear one elsewhere (best of the 24
colour permutations), then everything is repaired by tabucol."""
import itertools
import numpy as np
from flat import *
from cosets import cls

par = lambda x: bin(x).count("1") & 1


def res_bits(A):
    A = np.asarray(A, dtype=np.int64) % 2
    return (A * (1 << np.arange(12))).sum(axis=1)


def proper_pairs(dirs):
    """all row spaces {a,b,a+b} with (a.r,b.r) != (0,0) for every residue r of dirs; returns list of (a,b)"""
    R = sorted(set(int(x) for x in res_bits(dirs)))
    bitsR = np.array([[(r >> k) & 1 for k in range(12)] for r in R], dtype=np.int64)
    A = np.array([[(a >> k) & 1 for k in range(12)] for a in range(4096)], dtype=np.int64)
    dots = (A @ bitsR.T) & 1                    # (4096, |R|)
    zero = dots == 0
    Zf = zero.astype(np.float32)
    good = (Zf @ Zf.T) == 0                      # (a,b) good iff no r with a.r = b.r = 0
    spaces = set()
    for a, b in zip(*np.nonzero(good)):
        if a < b:
            spaces.add(frozenset((int(a), int(b), int(a) ^ int(b))))
    out = []
    for s in spaces:
        a, b = sorted(s)[:2]
        out.append((a, b))
    return sorted(out)


def lin_colour(pib, a, b):
    """pib: array of 12-bit residues of the points"""
    pa = np.array([par(int(x) & a) for x in pib]); pb = np.array([par(int(x) & b) for x in pib])
    return 2 * pa + pb


def fit_triangle(col, tri):
    """permute colours so that the triangle gets 0,1,2 (if its colours are distinct), else None"""
    c = [int(col[v]) for v in tri]
    if len(set(c)) < 3:
        return None
    rest = [x for x in range(4) if x not in c][0]
    perm = {c[0]: 0, c[1]: 1, c[2]: 2, rest: 3}
    return np.array([perm[int(x)] for x in col])


def conflicts(col, E):
    return int(np.sum(col[E[:, 0]] == col[E[:, 1]]))


def starts(P, E, tri, DP, CP, nbest=6):
    """candidate initial colourings sorted by number of conflicts"""
    pib = res_bits(P)
    C = cls(P)
    rowb = C[tri[0]][1]
    inrow = np.array([c[1] == rowb for c in C])
    cand = []
    linD = {ab: lin_colour(pib, *ab) for ab in DP}
    linC = {ab: lin_colour(pib, *ab) for ab in CP}
    for ab, col in list(linD.items()) + list(linC.items()):
        f = fit_triangle(col, tri)
        if f is not None:
            cand.append((conflicts(f, E), "lin", f))
    # mixed: row part D-linear, the rest conj-linear under the best permutation
    for ab, cD in linD.items():
        fD = fit_triangle(cD, tri)
        if fD is None:
            continue
        for ab2, cC in linC.items():
            best = None
            for perm in itertools.permutations(range(4)):
                pc = np.array(perm)[cC]
                col = np.where(inrow, fD, pc)
                k = conflicts(col, E)
                if best is None or k < best[0]:
                    best = (k, "mix", col)
            cand.append(best)
    cand.sort(key=lambda t: t[0])
    return cand[:nbest]


if __name__ == "__main__":
    import sys, json, time
    from tabu import tabucol
    D, CD, U = directions()
    DP, CP = proper_pairs(D), proper_pairs(CD)
    print("D-proper row spaces:", len(DP), " conj-proper:", len(CP), " common:", len(set(DP) & set(CP)))
    P = np.load(sys.argv[1]); st = json.load(open(sys.argv[2]))
    tri = st["triangle"]
    E, J = build_edges(P, U)
    t = time.time()
    S = starts(P, E, tri, DP, CP, nbest=5)
    print(f"n={len(P)} m={len(E)}; best starts (conflicts): {[(k, s) for k, s, _ in S]}  [{time.time() - t:.1f}s]")
    for k, s, col in S[:3]:
        t = time.time()
        c = tabucol(len(P), E, init=col, fixed=[(v, i) for i, v in enumerate(tri)], maxiter=3_000_000, seed=1)
        print(f"  start {s} with {k} conflicts -> tabucol {'OK' if c is not None else 'FAIL'} in {time.time() - t:.1f}s")
