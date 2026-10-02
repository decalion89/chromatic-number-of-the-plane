"""scan5.py d,d,... DMAX PMAX: for every D <= DMAX whose prime factors are <= PMAX, the unit vectors U_D of
Q(sqrt d)^2, the gates of gates2.py, and whether U_D closes a 5-cycle (five unit vectors with sum 0, so that the
unit-distance graph has cycles of length 5; there are no triangles). Prints each (d, D) with all gates open and a
5-cycle, then a summary per d. On 1-2 October the two quick successes of the field scan (d = 251 with D = 390,
d = 455 with D = 780) were exactly the prescreened pairs with a 5-cycle; it then chose the denominators for
d = 263, 299, 407, 599, 935 and 959 (263, 599 and 935 with D = 1020 and 959 with D = 2460 succeeded; 299 and
407 stalled)."""
import sys
import numpy as np
from units_fast import units_fast
from gates2 import report

OFF = 1 << 15


def keys(A):
    B = (A + OFF).astype(np.uint64)
    return (B[:, 0] << np.uint64(48)) | (B[:, 1] << np.uint64(32)) | (B[:, 2] << np.uint64(16)) | B[:, 3]


def has5(U):
    """some u1 + ... + u5 = 0 with every ui in U: a sum of three equals a sum of two (U = -U)"""
    S2 = np.unique((U[:, None, :] + U[None, :, :]).reshape(-1, 4), axis=0)
    ks = np.unique(keys(S2)); step = max(1, 2_000_000 // len(U))
    for i in range(0, len(S2), step):
        T = (S2[i:i + step][:, None, :] + U[None, :, :]).reshape(-1, 4)
        if np.isin(keys(T), ks).any():
            return True
    return False


def smooth(n, pmax):
    for p in range(2, pmax + 1):
        while n % p == 0:
            n //= p
    return n == 1


if __name__ == "__main__":
    ds = [int(x) for x in sys.argv[1].split(",")]; DMAX, PMAX = int(sys.argv[2]), int(sys.argv[3])
    for d in ds:
        hits = []
        for D in range(2, DMAX + 1):
            if not smooth(D, PMAX):
                continue
            U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
            if len(U) < 40 or len(U) > 400:
                continue
            r = report(d, D); closed = [g for g in ("2", "3", "5") if r.get(g) == 0]
            if closed:
                continue
            if has5(U):
                hits.append((D, len(U)))
                print(f"d={d} D={D}: {len(U)} directions, gates {r}: closed walk of length 5", flush=True)
        print(f"SUMMARY d={d}: {len(hits)} denominators with a 5-walk and all gates open: {hits[:30]}", flush=True)
