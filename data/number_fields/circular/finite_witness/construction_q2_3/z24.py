"""Exact arithmetic in Q(zeta_24) on integer coordinate vectors over the power basis 1, zeta, ..., zeta^7
(zeta^8 = zeta^4 - 1), with complex conjugation zeta -> zeta^-1; a point of Q(sqrt2, sqrt3)^2 is x + iy."""
import numpy as np


def pmul(a, b):
    r = [0] * 15
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    r[i + j] += x * y
    for k in range(14, 7, -1):
        c = r[k]
        if c:
            r[k] = 0; r[k - 4] += c; r[k - 8] -= c
    return r[:8]


def zpow(j):
    e = [1, 0, 0, 0, 0, 0, 0, 0]
    for _ in range(j % 24):
        e = pmul(e, [0, 1, 0, 0, 0, 0, 0, 0])
    return e


CONJ = [zpow(-k) for k in range(8)]            # coordinates of conj(zeta^k)


def conj(a):
    r = [0] * 8
    for k, x in enumerate(a):
        if x:
            for t in range(8):
                r[t] += x * CONJ[k][t]
    return r


def norm2(a):
    """a * conj(a) as a coordinate vector (lies in the real subfield)"""
    return pmul(a, conj(a))


def is_unit(diff, D):
    return norm2(list(diff)) == [D * D, 0, 0, 0, 0, 0, 0, 0]


ZC = np.exp(2j * np.pi * np.arange(8) / 24)


def complex_of(P, D):
    return (np.asarray(P, dtype=float) @ ZC) / D
