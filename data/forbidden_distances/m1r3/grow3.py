"""Greedy growth of a 3D configuration for the inclusion-exclusion program (floating point exploration).
usage: python3 grow3.py SEED.npy OUT_PREFIX NMAX [MAXEVAL] [MINDEG]"""
import sys, itertools, time, json
import numpy as np, os
from multiprocessing import Pool
NPROC = int(os.environ.get('NPROC', 2))
from ielp3 import solve, unit_graph, indep_sets
cuts = np.arange(0.10, 0.20, 0.0025)
X = np.load(sys.argv[1]); pre = sys.argv[2]; nmax = int(sys.argv[3])
maxeval = int(sys.argv[4]) if len(sys.argv) > 4 else 60; mindeg = int(sys.argv[5]) if len(sys.argv) > 5 and not sys.argv[5].startswith('--') else 3
def sphere3(a, b, c):
    """points at distance 1 from a, b, c"""
    ex = b - a; d = np.linalg.norm(ex)
    if d < 1e-9 or d > 2: return []
    ex /= d; i = ex @ (c - a); ey = c - a - i * ex; ny = np.linalg.norm(ey)
    if ny < 1e-9: return []
    ey /= ny; ez = np.cross(ex, ey); j = ey @ (c - a)
    x = d / 2; y = (i * i + j * j - 2 * i * x) / (2 * j)
    z2 = 1 - x * x - y * y
    if z2 < -1e-12: return []
    z = np.sqrt(max(z2, 0))
    return [a + x * ex + y * ey + s * z * ez for s in ((1, -1) if z > 1e-9 else (1,))]
RAT = '--rational' in sys.argv
from fractions import Fraction as Fr
def rational(p, den=2000):
    q = p * np.sqrt(2)
    return all(abs(float(Fr(float(x)).limit_denominator(den)) - x) < 1e-11 for x in q)
def value(X):
    D, adj = unit_graph(X); M = indep_sets(len(X), adj)
    v, dt, info = solve(X, M=M, adj=adj, cuts=cuts); return v, len(M), dt
v0, na, dt = value(X); print(f"seed n={len(X)} value {v0:.6f} atoms {na} ({dt:.1f}s)", flush=True)
while len(X) < nmax:
    T = time.time(); cand = []
    for a, b, c in itertools.combinations(range(len(X)), 3):
        for p in sphere3(X[a], X[b], X[c]):
            dist = np.linalg.norm(X - p, axis=1)
            if dist.min() < 1e-6: continue
            if RAT and not rational(p): continue
            cand.append(p)
    if not cand: print('no candidates'); break
    C = np.array(cand); key = np.round(C, 6) + 0.0
    _, idx = np.unique(key, axis=0, return_index=True); C = C[np.sort(idx)]
    deg = np.array([(np.abs(np.linalg.norm(X - p, axis=1) - 1) < 1e-9).sum() for p in C])
    order = np.lexsort((np.linalg.norm(C, axis=1), -deg)); order = [o for o in order if deg[o] >= mindeg][:maxeval]
    print(f"n={len(X)}: {len(C)} candidates, degrees {np.bincount(deg).tolist()}, evaluating {len(order)}", flush=True)
    with Pool(NPROC) as pool:
        res = pool.map(value, [np.vstack([X, C[o]]) for o in order])
    k = int(np.argmin([r[0] for r in res])); v, na = res[k][0], res[k][1]; o = order[k]
    top = sorted(range(len(res)), key=lambda i: res[i][0])[:5]
    print('   top5', [(round(res[i][0], 6), int(deg[order[i]])) for i in top], flush=True)
    X = np.vstack([X, C[o]]); np.save(f"{pre}_{len(X)}.npy", X)
    print(f"  added {np.round(C[o], 5).tolist()} deg {deg[o]}: value {v:.6f} atoms {na} ({time.time() - T:.0f}s)", flush=True)
    if na > 4e5: break
