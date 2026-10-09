"""For each d and each D in a list: the number of directions, the gates of gates2.py (a closed gate colours every graph
built from U_D), and the shortest odd closed walk (up to 7). One line per (d, D); the last line per d names the
best D: all gates open, shortest odd walk, then most directions (at most NMAX). usage: prescreen.py d,d,... D,D,..."""
import sys, time
import numpy as np
from units_fast import units_fast
from gates2 import report
OFF = 1 << 15
def keys(A):
    B = (A + OFF).astype(np.uint64)
    return (B[:, 0] << np.uint64(48)) | (B[:, 1] << np.uint64(32)) | (B[:, 2] << np.uint64(16)) | B[:, 3]
def shortest_odd(U, kmax=3):
    S = np.zeros((1, 4), dtype=np.int64)
    for k in range(kmax + 1):
        ks = np.unique(keys(S)); step = max(1, 2_000_000 // len(U))
        for i in range(0, len(S), step):
            T = (S[i:i + step][:, None, :] + U[None, :, :]).reshape(-1, 4)
            if np.isin(keys(T), ks).any():
                return 2 * k + 1
        if k < kmax:
            S = (S[:, None, :] + U[None, :, :]).reshape(-1, 4)
            _, idx = np.unique(keys(S), return_index=True); S = S[idx]
    return None
ds = [int(x) for x in sys.argv[1].split(",")]; Ds = [int(x) for x in sys.argv[2].split(",")]
NMAX = 260
for d in ds:
    best = []
    for D in Ds:
        U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
        if len(U) < 60:
            continue
        r = report(d, D); closed = [g for g in ("2", "3", "5") if r.get(g) == 0]
        L = shortest_odd(U) if not closed and len(U) <= NMAX else None
        print(f"d={d} D={D}: {len(U)} directions, gates {r}{' CLOSED ' + ','.join(closed) if closed else ''}, shortest odd {L}", flush=True)
        if not closed and L is not None and len(U) <= NMAX:
            best.append((L, -len(U), D))
    best.sort()
    print(f"BEST d={d}: {[(D, -n, L) for L, n, D in best[:4]]}", flush=True)
