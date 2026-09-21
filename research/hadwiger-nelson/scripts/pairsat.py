"""Every pair at one given squared distance, by the same bucketing the unit
edge finder uses.

The scans in this work computed the pairs at a second distance with a Python
loop over vertices, one vectorised row-difference each.  That is fine at 1581
points and hopeless at nineteen thousand: the loop overhead dominates and the
work is quadratic on top of it.

The unit-distance finder does not do that.  It buckets the points into cells
of side one, so that a pair at distance one lies in the same cell or in one of
the eight around it, then confirms every float candidate by exact integer
field arithmetic.  The same argument holds for any distance with cells of that
side, so the same routine serves -- with the cell side and the target changed.
"""
import numpy as np


def pairs_at(basis, rows, target_num, eps=1e-6):
    """Exact pairs whose squared distance equals target_num / basis.D**2."""
    n = len(rows)
    if n == 0:
        return []
    xy = basis.floats(rows)
    D2 = basis.D * basis.D
    side = (float(target_num) / D2) ** .5
    if side <= 0:
        return []
    cell = np.floor(xy / side).astype(np.int64)
    span = int(cell.max() - cell.min()) + 3
    base = cell.min()
    key = (cell[:, 0] - base) * span + (cell[:, 1] - base)
    order = np.argsort(key, kind="stable")
    ks = key[order]
    uniq, starts = np.unique(ks, return_index=True)
    ends = np.append(starts[1:], len(ks))
    block = {int(k): (int(s), int(e)) for k, s, e in zip(uniq, starts, ends)}
    dim = basis.dim
    out = []
    tgt = float(target_num) / D2
    for k, (s, e) in block.items():
        ia = order[s:e]
        pa = xy[ia]
        for dx, dy in ((0, 0), (1, 0), (0, 1), (1, 1), (-1, 1)):
            k2 = k + dx * span + dy
            if k2 not in block:
                continue
            s2, e2 = block[k2]
            ib = order[s2:e2]
            pb = xy[ib]
            d = pa[:, None, :] - pb[None, :, :]
            dd = d[:, :, 0] ** 2 + d[:, :, 1] ** 2
            m = np.abs(dd - tgt) < eps * max(1.0, tgt)
            if dx == 0 and dy == 0:
                m &= ia[:, None] < ib[None, :]
            ii, jj = np.nonzero(m)
            if not len(ii):
                continue
            A, B = ia[ii], ib[jj]
            diff = rows[A] - rows[B]
            sq = (basis._field_square(diff[:, :dim])
                  + basis._field_square(diff[:, dim:]))
            good = sq[:, 0] == target_num
            for c in range(1, dim):
                good &= sq[:, c] == 0
            for a, b in zip(A[good], B[good]):
                out.append((int(min(a, b)), int(max(a, b))))
    return sorted(set(out))
