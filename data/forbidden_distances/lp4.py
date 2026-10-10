"""Clique-strengthened LP with local refinement of the placements (float screening; certify with certify_lp2.py).
usage: python3 lp4.py "d1;d2" shape1,shape2,... out.json
Placements: origin p = centroid + (ox, oy), first on the grid (1/2)Z^2 within radius 5/2, then refined with step 1/10
around the active ones (two rounds). Offsets are stored as exact fractions."""
import sys, json, itertools, numpy as np
from fractions import Fraction as Fr
from scipy.special import j0
from scipy.optimize import linprog
from shapes import SHAPES, ev
RMIN = 0.1

def build(D, shapes, placements, t):
    cols = [np.ones_like(t)] + [j0(d * t) for d in D]
    m0 = [1.0] * (1 + len(D)); obj = [1.0] + [0.0] * len(D); bounds = [(0, None)] + [(None, None)] * len(D)
    for (name, ox, oy) in placements:
        W = shapes[name]; cen = W.mean(axis=0)
        p = cen + np.array([float(ox), float(oy)])
        r = np.linalg.norm(W - p, axis=1)
        cols.append(sum(j0(ri * t) for ri in r)); m0.append(float(len(W))); obj.append(1.0); bounds.append((0, None))
    A = np.column_stack(cols)
    A_ub = np.vstack([-A, -np.array(m0)[None, :]]); b_ub = np.concatenate([np.zeros(len(t)), [-1.0]])
    return linprog(np.array(obj), A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

def ok(shapes, name, ox, oy):
    W = shapes[name]; p = W.mean(axis=0) + np.array([float(ox), float(oy)])
    return np.linalg.norm(W - p, axis=1).min() > RMIN

def main():
    D = [ev(e) for e in sys.argv[1].split(';')]
    names = sys.argv[2].split(','); out = sys.argv[3]
    shapes = {n: np.array([[ev(x), ev(y)] for x, y in SHAPES[n]]) for n in names}
    for n in names:
        W = shapes[n]
        for a, b in itertools.combinations(range(len(W)), 2):
            dd = np.linalg.norm(W[a] - W[b]); assert min(abs(dd - d) for d in D) < 1e-9, (n, dd)
    t = np.concatenate([np.linspace(1e-6, 30.0, 12000), np.linspace(30.0, 400.0, 30000)])
    g = Fr(1, 2)
    pl = [(n, i * g, j * g) for n in names for i in range(-5, 6) for j in range(-5, 6) if ok(shapes, n, i * g, j * g)]
    for rnd, step in enumerate([Fr(1, 10), Fr(1, 50)]):
        res = build(D, shapes, pl, t)
        x = res.x; k0 = 1 + len(D)
        act = [pl[k] for k in range(len(pl)) if x[k0 + k] > 1e-9]
        print(f'round {rnd}: {len(pl)} placements, LP {res.fun:.6f} (1/LP {1/res.fun:.4f}), active {len(act)}', flush=True)
        new = set(act)
        for (n, ox, oy) in act:
            for a in range(-5, 6):
                for b in range(-5, 6):
                    q = (n, ox + a * step, oy + b * step)
                    if ok(shapes, *q):
                        new.add(q)
        pl = sorted(new)
    res = build(D, shapes, pl, t); x = res.x; k0 = 1 + len(D)
    act = [{'shape': pl[k][0], 'ox': str(pl[k][1]), 'oy': str(pl[k][2]), 'y': float(x[k0 + k])}
           for k in range(len(pl)) if x[k0 + k] > 1e-9]
    print(f'final: {len(pl)} placements, LP {res.fun:.6f} (1/LP {1/res.fun:.4f}), active {len(act)}', flush=True)
    json.dump({'D': sys.argv[1], 'z': [float(v) for v in x[:k0]], 'cliques': act, 'lp': float(res.fun)},
              open(out, 'w'), indent=1)

if __name__ == '__main__':
    main()
