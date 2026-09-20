"""Forced-pair search that does not enumerate the pairs first.

At 27000 points there are 360 million candidate pairs, and building that list
costs more than the search.  But the filter throws almost all of them away
immediately, so the list never needs to exist: sample the colourings FIRST,
then walk the pairs once and keep only those that agree in every sample.  Ten
samples leave about a ten-millionth of them, so the survivors fit in a
handful of tuples and memory stays flat.

Vectorised, the agreement test is a single numpy reduction per row, and the
geometry is only computed for the rows that survive it -- which is almost
none.  What would have been ten minutes of Python is a few seconds.
"""
import time
import numpy as np
from fractions import Fraction as Fr
from pysat.solvers import Solver
from hn.homcol import closable_distance


def sample_colourings(cls, n, k, samples=14, seed=0):
    """Distinct proper k-colourings, via randomised saved phases.

    Cadical ignores set_phases and returns one colouring however often it is
    asked; glucose honours them.  Diversity is the whole point, so glucose.
    """
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None
    rng = np.random.default_rng(seed)
    out = np.empty((samples, n), dtype=np.int8)
    for s in range(samples):
        sv.set_phases([int(v) for v in
                       (rng.integers(0, 2, n * k) * 2 - 1)
                       * np.arange(1, n * k + 1)])
        sv.solve()
        m = sv.get_model()
        arr = np.frombuffer(np.array(m[:n * k], dtype=np.int64), dtype=np.int64)
        out[s] = (arr.reshape(n, k) > 0).argmax(axis=1)
    sv.delete()
    return out


def survivors(cols, zf, denom=1584, lim=36.0, report=None):
    """Pairs that agree in every sampled colouring and sit at a closable
    rational distance.  Geometry is only touched for pairs that agree."""
    n = cols.shape[1]
    xs = np.array([p[0] for p in zf])
    ys = np.array([p[1] for p in zf])
    out, t0 = [], time.time()
    for i in range(n - 1):
        agree = (cols[:, i + 1:] == cols[:, i:i + 1]).all(axis=0)
        idx = np.nonzero(agree)[0]
        if idx.size:
            j = idx + i + 1
            v = (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2
            for jj, vv in zip(j, v):
                if vv > lim:
                    continue
                D = Fr(round(vv * denom), denom)
                if abs(float(D) - vv) > 1e-7 or D == 1:
                    continue
                if closable_distance(D):
                    out.append((i, int(jj)))
        if report and i % report == 0:
            print(f"    ... {i}/{n}, {len(out)} survivors "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    return out


def forced_among(cls, k, cands):
    sv = Solver(name="cd19", bootstrap_with=cls)
    hits = [(i, j) for i, j in cands
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    sv.delete()
    return hits
