"""Greedy extension of a point set for the one-distance inclusion-exclusion program (floating point, exploration).
Candidates: points at distance 1 from at least two points of X and at distance >= 0.1 from all of them (the rule of
Ambrus et al.). For each candidate c the program for X + {c} is solved; results go to a JSON-lines file.
usage: python3 explore.py BASE.json|BASE.npy OUT.jsonl WORKER NWORKERS [tmax step mindegree]"""
import sys, json, time
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ielp

def canon_fast(P, nd=6):
    """canonical form under isometries, using only frames on diameter pairs (the diameter is an isometry invariant)"""
    k = len(P)
    D = np.abs(P[:, None] - P[None, :]); dmax = D.max()
    frames = [(a, b) for a in range(k) for b in range(k) if a != b and D[a, b] > dmax - 1e-7]
    best = None
    for a, b in frames:
        v = P[b] - P[a]; w = (P - P[a]) * np.conj(v) / abs(v)
        for z in (w, np.conj(w)):
            key = tuple(sorted((round(x.real, nd) + 0.0, round(x.imag, nd) + 0.0) for x in z))
            if best is None or key < best: best = key
    return best

def congruences_fast(X, M):
    classes = {}
    for m in M:
        m = int(m)
        if m & (m - 1) == 0: continue
        idx = [i for i in range(len(X)) if m >> i & 1]
        classes.setdefault(canon_fast(X[idx]), []).append(m)
    rel = [(ms[0], m) for ms in classes.values() for m in ms[1:]]
    return rel, len(classes)

def candidates(X, mind=0.1):
    out = {}
    n = len(X)
    for i in range(n):
        for j in range(i + 1, n):
            d = X[j] - X[i]; r = abs(d)
            if r >= 2 - 1e-12 or r < 1e-12: continue
            h = np.sqrt(1 - r * r / 4); m = (X[i] + X[j]) / 2
            for s in (1, -1):
                c = m + s * 1j * d / r * h
                if np.min(np.abs(X - c)) < mind: continue
                key = (round(c.real, 8), round(c.imag, 8))
                out.setdefault(key, c)
    cs = list(out.values())
    deg = [int(np.sum(np.abs(np.abs(X - c) - 1) < 1e-9)) for c in cs]
    return cs, deg

if __name__ == '__main__':
    if sys.argv[1].endswith('.json'):                       # an exact configuration, as in configurations/
        import exactconf as ec
        pts, _ = ec.load(sys.argv[1])
        base = np.array([complex(float(ec.kval(re)), float(ec.kval(im))) for re, im in pts])
    else:
        base = np.load(sys.argv[1])
    outp = sys.argv[2]; w, nw = int(sys.argv[3]), int(sys.argv[4])
    tmax = float(sys.argv[5]) if len(sys.argv) > 5 else 100.0; step = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
    dmin = int(sys.argv[7]) if len(sys.argv) > 7 else 2
    cs, deg = candidates(base)
    order = sorted((k for k in range(len(cs)) if deg[k] >= dmin), key=lambda k: -deg[k])   # many unit neighbours first
    mine = [k for r, k in enumerate(order) if r % nw == w]
    print(f'{len(cs)} candidates, {len(order)} of degree >= {dmin}; worker {w}/{nw} takes {len(mine)}', flush=True)
    for k in mine:
        c = cs[k]; X = np.append(base, c)
        t0 = time.time()
        D, adj = ielp.unit_graph(X); M = ielp.indep_sets(len(X), adj); rel, ncl = congruences_fast(X, M)
        v, dt, info = ielp.solve(X, tmax, step, M=M, rel=rel)
        rec = dict(k=k, re=c.real, im=c.imag, deg=deg[k], value=v, atoms=len(M), rel=len(rel), lp_s=round(dt, 1),
                   total_s=round(time.time() - t0, 1), status=info['status'])
        with open(outp, 'a') as f: f.write(json.dumps(rec) + '\n')
        print(json.dumps(rec), flush=True)
