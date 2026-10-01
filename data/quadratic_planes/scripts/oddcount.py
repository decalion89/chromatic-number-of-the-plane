"""How many points lie both at 3 steps and at 4 steps from the origin (each gives odd closed walks of length 7),
and at 2 and 3 steps (length 5), in Cay(Z^4, U_D) over Q(sqrt d)? usage: oddcount.py d D"""
import sys
import numpy as np
from units_fast import units_fast
d, D = int(sys.argv[1]), int(sys.argv[2])
U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
OFF = 1 << 15
def keys(A):
    B = (A + OFF).astype(np.uint64)
    return (B[:, 0] << np.uint64(48)) | (B[:, 1] << np.uint64(32)) | (B[:, 2] << np.uint64(16)) | B[:, 3]
def step(S):
    T = (S[:, None, :] + U[None, :, :]).reshape(-1, 4)
    k, idx = np.unique(keys(T), return_index=True)
    return T[idx], k
S0 = np.zeros((1, 4), dtype=np.int64)
S1, k1 = step(S0); S2, k2 = step(S1); S3, k3 = step(S2)
c5 = np.intersect1d(k2, k3).size
c7 = 0
step_n = max(1, 2_000_000 // len(U))
hits = []
for i in range(0, len(S3), step_n):
    T = (S3[i:i + step_n][:, None, :] + U[None, :, :]).reshape(-1, 4)
    kt = np.unique(keys(T))
    hits.append(kt[np.isin(kt, k3)])
c7 = np.unique(np.concatenate(hits)).size if hits else 0
print(f"d={d} D={D}: {len(U)} directions; |S2| = {len(k2)}, |S3| = {len(k3)}; points at 2 and 3 steps: {c5}; at 3 and 4 steps: {c7}", flush=True)
