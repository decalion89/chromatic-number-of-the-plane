"""limit_test.py: the N -> infinity limit test for a finite set V of unit vectors of F^2 (F any number field, nf.Field).
By Proposition 1 and Lemma 5 of notes/four_colours_11_mod_12.md, for large N = 5^k the set G_N * V has a character into
[1/3, 2/3] only if some Q(i)-linear phi : L -> C has phi(v) in E for all v in V, E = ((1+i)/2 + Z[i]) u
((a+bi)/3 + Z[i], a, b in {1, 2}).  With z = 6 phi(v) in Z^2: type c <=> z = (1,1) mod 2 and (0,0) mod 3;
type q <=> z = (0,0) mod 2 and both coordinates nonzero mod 3.  L = (image of phi) cap Z^{2m}; we test L/2L x L/3L."""
import sys, itertools
from fractions import Fraction as Fr
from math import lcm
import numpy as np


def rows_of(Fd, V):
    M = []
    for (X, Y) in V:
        xv, yv = Fd.vec(X), Fd.vec(Y)
        re, im = [], []
        for b in range(Fd.n):
            re += [xv[b], -yv[b]]
            im += [yv[b], xv[b]]
        M.append(re); M.append(im)
    return M


def lattice(Fd, V):
    """basis of (Q-span of the image of phi) cap Z^{2m}, exactly, with PARI/GP matrixqz(A, -2)"""
    import subprocess, tempfile, os
    M = rows_of(Fd, V)
    D = 1
    for r in M:
        for x in r: D = lcm(D, Fr(x).denominator)
    A = [[int(Fr(x) * D) for x in r] for r in M]
    mat = "[" + ";".join(",".join(str(x) for x in row) for row in A) + "]"
    script = "A = " + mat + ";\nB = matrixqz(A, -2);\nprint(#B);\nfor(j=1,#B, print(Vec(B[,j])));\nquit;\n"
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as fh:
        fh.write(script)
        fn = fh.name
    out = subprocess.run(["gp", "-q", "-s", "400000000", fn], capture_output=True, text=True, check=True).stdout.split("\n")
    os.unlink(fn)
    r = int(out[0]); L = []
    for j in range(r):
        L.append([int(x) for x in out[1 + j].strip()[1:-1].split(",")])
    return L


def feasible(Fd, V, want=False):
    L = lattice(Fd, V)
    r = len(L); m = len(V)
    B2m = np.array([[x % 2 for x in row] for row in L], dtype=np.int64)
    B3m = np.array([[x % 3 for x in row] for row in L], dtype=np.int64)
    # mod 2
    c2 = np.array(list(itertools.product(range(2), repeat=r)), dtype=np.int64)
    X2 = (c2 @ B2m) % 2
    pairs = X2.reshape(len(c2), m, 2)
    is11 = (pairs[:, :, 0] == 1) & (pairs[:, :, 1] == 1)
    is00 = (pairs[:, :, 0] == 0) & (pairs[:, :, 1] == 0)
    okx = np.all(is11 | is00, axis=1)
    if not okx.any():
        return (False, None) if want else False
    pats = []
    seen = set()
    for k in np.nonzero(okx)[0]:
        pat = tuple(is11[k])
        if pat not in seen:
            seen.add(pat); pats.append(np.array(pat))
    # enumerate L/3L in chunks (memory), test every surviving mod-2 pattern
    import itertools as it
    allc = it.product(range(3), repeat=r)
    while True:
        chunk = list(it.islice(allc, 40000))
        if not chunk: break
        c3 = np.array(chunk, dtype=np.int64)
        Y3 = (c3 @ B3m) % 3
        yp = Y3.reshape(len(c3), m, 2)
        y00 = (yp[:, :, 0] == 0) & (yp[:, :, 1] == 0)
        ynz = (yp[:, :, 0] != 0) & (yp[:, :, 1] != 0)
        for cmask in pats:
            good = np.all(np.where(cmask[None, :], y00, ynz), axis=1)
            if good.any():
                return (True, "".join("c" if t else "q" for t in cmask)) if want else True
    return (False, None) if want else False
