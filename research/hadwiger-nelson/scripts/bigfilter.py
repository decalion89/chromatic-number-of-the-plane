"""Forced-pair search that does not enumerate the pairs first.

At 27000 points there are 360 million candidate pairs, and building that list
costs more than the search.  But the filter throws almost all of them away
immediately, so the list never needs to exist: sample the colourings FIRST,
then walk the pairs once and keep only those that agree in every sample.  Ten
samples leave about a ten-millionth of them, so the survivors fit in a
handful of tuples and memory stays flat.

Even that is too slow at scale: a reduction per row is O(samples * n^2), ten
billion element comparisons at 27000 points, three hours a union.  But two
vertices agree in EVERY sample exactly when their colour signatures across the
samples are identical, so the pairs never have to be compared at all -- bucket
the vertices by signature and the candidates are the pairs inside a bucket.
Fourteen samples give 5^14 signatures, six billion of them, so almost every
bucket is a singleton and the few collisions are the entire candidate set.
O(n) instead of O(n^2), and the size of the graph stops mattering.
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
    rational distance.

    Vertices are bucketed by their colour signature across the samples, which
    is what "agrees in every sample" means, so no pair is ever compared.  The
    geometry is then computed only inside buckets, of which almost all are
    singletons.
    """
    n = cols.shape[1]
    key = np.zeros(n, dtype=np.int64)
    for s in range(cols.shape[0]):
        key = key * 5 + cols[s].astype(np.int64)
    order = np.argsort(key, kind="stable")
    out, start = [], 0
    ks = key[order]
    while start < n:
        stop = start + 1
        while stop < n and ks[stop] == ks[start]:
            stop += 1
        if stop - start > 1:
            grp = order[start:stop]
            for a in range(len(grp)):
                for b in range(a + 1, len(grp)):
                    i, j = int(grp[a]), int(grp[b])
                    v = ((zf[i][0] - zf[j][0]) ** 2
                         + (zf[i][1] - zf[j][1]) ** 2)
                    if v > lim:
                        continue
                    D = Fr(round(v * denom), denom)
                    if abs(float(D) - v) > 1e-7 or D == 1:
                        continue
                    if closable_distance(D):
                        out.append((min(i, j), max(i, j)))
        start = stop
    if report:
        print(f"    {n} vertices, {n - len(np.unique(key))} in collisions, "
              f"{len(out)} survivors", flush=True)
    return sorted(set(out))


def forced_among(cls, k, cands):
    sv = Solver(name="cd19", bootstrap_with=cls)
    hits = [(i, j) for i, j in cands
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    sv.delete()
    return hits
