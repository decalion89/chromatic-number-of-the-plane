"""Exact shortest odd closed walk in Cay(Z^4, U_D) over Q(sqrt d), by meeting in the middle: the shortest odd
closed walk has length 2k+1 for the least k such that some sum of k directions plus one more direction is again
a sum of k directions (S_k + U meets S_k). usage: oddmitm.py d D [kmax]"""
import sys, time
import numpy as np
from units_fast import units_fast
d, D = int(sys.argv[1]), int(sys.argv[2]); kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 4
U = np.asarray(units_fast(d, D), dtype=np.int64).reshape(-1, 4)
OFF = 1 << 15
def keys(A):
    assert np.all(np.abs(A) < OFF - 1)
    B = (A + OFF).astype(np.uint64)
    return (B[:, 0] << np.uint64(48)) | (B[:, 1] << np.uint64(32)) | (B[:, 2] << np.uint64(16)) | B[:, 3]
t = time.time()
S = np.zeros((1, 4), dtype=np.int64)                      # S_0
res = None
for k in range(0, kmax + 1):
    ks = np.unique(keys(S))
    hit = False
    for i in range(0, len(S), max(1, 2_000_000 // len(U))):
        T = (S[i:i + max(1, 2_000_000 // len(U))][:, None, :] + U[None, :, :]).reshape(-1, 4)
        if np.isin(keys(T), ks, assume_unique=False).any():
            hit = True; break
    if hit:
        res = 2 * k + 1; break
    if k == kmax:
        break
    S = (S[:, None, :] + U[None, :, :]).reshape(-1, 4)
    _, idx = np.unique(keys(S), return_index=True); S = S[idx]
print(f"d={d} D={D}: {len(U)} directions; shortest odd closed walk: {res if res else f'> {2 * kmax + 1}'}  [{time.time() - t:.0f}s, |S_k| = {len(S)}]", flush=True)
