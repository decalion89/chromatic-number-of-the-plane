"""The three-point bound for the finite planes (notes/local_colourings.md, section 14; scripts/threepoint.py).

Structural tests, with no SDP solver:
- the rotation blocks are the compressions of the explicit matrices M1 and M0 by an orthonormal basis
  (so the blocks are psd exactly when the matrices are), and the localizing blocks are compressions of the
  explicit localizing matrices;
- an actual independent set satisfies every constraint, with objective |S|.
The stored certificates are checked in tests/test_threepoint_certificates.py (marked slow).
"""
import os, sys, json
import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.append(os.path.join(HERE, '..', 'scripts'))   # appended: scripts/six.py must not shadow the six package

from threepoint import build          # noqa: E402


def explicit_matrices(D, v):
    q = D['meta']['q']; n = q * q
    cls, lab, zid, unit, gcol = D['cls'], D['lab'], D['zid'], D['unit'], D['gcol']
    X, Y = np.divmod(np.arange(n), q)
    idx = lambda x, y: (x % q) * q + (y % q)
    g = lambda w: 1.0 if w == 0 else (0.0 if cls[w] == 1 else v[gcol[cls[w]]])

    def z(a, b):
        if a == b: return g(a)
        if a == 0: return g(b)
        if b == 0: return g(a)
        p = a * n + b
        return 0.0 if unit[p] else v[zid[lab[p]]]
    sub = lambda u, w: int(idx(X[u] - X[w], Y[u] - Y[w]))
    M1 = np.array([[z(a, b) for b in range(n)] for a in range(n)])
    M0 = np.array([[(1 - g(a)) if a == b else g(sub(b, a)) - z(a, b) for b in range(1, n)] for a in range(1, n)])
    return M1, M0, g


def block_values(ent, size, v):
    B = np.zeros((size, size))
    for (i, j), (acc, cst) in ent.items():
        B[i, j] = B[j, i] = np.real(cst) + sum(np.real(c) * v[vi] for vi, c in acc.items())
    return B


@pytest.mark.parametrize("q,kind", [(7, 'std'), (11, 'std'), (13, 'std'), (13, 'inert')])
def test_rotation_blocks_are_orthonormal_compressions(q, kind):
    n = q * q; K0 = n // 4
    D = build(q, kind, verbose=False, K0=K0)
    v = np.random.default_rng(q).random(D['nv'])
    M1, M0, g = explicit_matrices(D, v)
    meta = D['meta']; N = meta['N']; m = meta['m']
    X, Y = np.divmod(np.arange(n), q)
    sig = np.array(meta['sig'])
    s_img = (sig[0, 0] * X + sig[0, 1] * Y) % q * q + (sig[1, 0] * X + sig[1, 1] * Y) % q
    powimg = np.zeros((N, n), int); powimg[0] = np.arange(n)
    for j in range(1, N): powimg[j] = s_img[powimg[j - 1]]
    reps, partner, shift = meta['reps'], meta['partner'], meta['shift']
    omega = np.exp(2j * np.pi / N)
    dims = {1: 0, 0: 0}; bi = 0
    for which in (1, 0):
        for k in range(N // 2 + 1):
            e = []
            for a in range(m):
                f = np.zeros(n, complex)
                for j in range(N): f[powimg[j, reps[a]]] += omega ** (k * (j - shift[a])) / np.sqrt(N)
                e.append(f)
            vecs = []; seen = set()
            for a in range(m):
                if a in seen: continue
                b = partner[a]
                if b < 0: vecs.append(e[a]); seen.add(a)
                else:
                    vecs += [(e[a] + e[b]) / np.sqrt(2), 1j * (e[a] - e[b]) / np.sqrt(2)]; seen.update((a, b))
            if which == 1 and k == 0:
                e0 = np.zeros(n, complex); e0[0] = 1; vecs = [e0] + vecs
            dims[which] += len(vecs) * (1 if (k == 0 or 2 * k == N) else 2)
            if which == 1:
                V, Mx = np.array(vecs), M1
            else:
                V, Mx = np.array(vecs)[:, 1:], M0
                if k == 0:
                    Mx = np.zeros((n, n)); Mx[1:, 1:] = M0; Mx[0, 0] = n / K0 - 1
                    for a in range(1, n): Mx[0, a] = Mx[a, 0] = 1 - g(a)
                    ee = np.zeros(n); ee[0] = 1
                    V = np.vstack([ee, np.hstack([np.zeros((len(vecs), 1)), V])])
            size, ent = D['blocks'][bi]; bi += 1
            B = V.conj() @ Mx @ V.T
            assert np.abs(B - block_values(ent, size, v)).max() < 1e-9 and np.abs(B.imag).max() < 1e-9
    assert dims == {1: n, 0: n - 1}


@pytest.mark.parametrize("q,kind,opt", [(13, 'std', 'trilocal'), (13, 'std', 'edge'), (11, 'std', 'edge')])
def test_localizing_blocks_are_compressions(q, kind, opt):
    n = q * q; K0 = n // 4
    D = build(q, kind, verbose=False, K0=K0, **{opt: True})
    vecsets = D['localizing_vectors']
    v = np.random.default_rng(1).random(D['nv'])
    M1, M0, g = explicit_matrices(D, v)
    cls, lab, zid, unit = D['cls'], D['lab'], D['zid'], D['unit']
    X, Y = np.divmod(np.arange(n), q)
    sub = lambda u, w: int(((X[u] - X[w]) % q) * q + (Y[u] - Y[w]) % q)
    U = np.nonzero(cls == 1)[0]; Us = set(U.tolist()); e = int(U[0])
    T = [0, e] if opt == 'edge' else [0, e, next(int(u) for u in U if sub(u, e) in Us)]

    def z3(t, a, b):
        p = sub(a, t) * n + sub(b, t)
        return 0.0 if unit[p] else v[zid[lab[p]]]
    pts = [p for p in range(n) if p not in T]
    ix = {p: i + 1 for i, p in enumerate(pts)}
    L = np.zeros((len(pts) + 1, len(pts) + 1)); L[0, 0] = n / K0 - len(T)
    for a in pts:
        L[0, ix[a]] = L[ix[a], 0] = L[ix[a], ix[a]] = 1 - sum(g(sub(a, t)) for t in T)
        for b in pts:
            if b != a: L[ix[a], ix[b]] = g(sub(b, a)) - sum(z3(t, a, b) for t in T)
    for (size, ent), vecs in zip(D['blocks'][D['meta']['nrot']:], vecsets):
        off = size - len(vecs)
        F = np.zeros((size, len(pts) + 1))
        if off: F[0, 0] = 1
        for i, coefs in enumerate(vecs):
            for p, c in coefs.items(): F[i + off, ix[p]] += c
        assert np.abs(F @ L @ F.T - block_values(ent, size, v)).max() < 1e-9


def large_independent_set(q, kind, tries=300, seed=0):
    """the largest of many random greedy independent sets containing 0 (any independent set will do)"""
    from threepoint import form
    n0 = form(q, kind)
    U = [(a, b) for a in range(q) for b in range(q) if (a * a + n0 * b * b) % q == 1]
    nb = [[((i // q + a) % q) * q + (i % q + b) % q for a, b in U] for i in range(q * q)]
    rng = np.random.default_rng(seed)
    best = []
    for _ in range(tries):
        order = [0] + list(rng.permutation(np.arange(1, q * q)))
        blocked = np.zeros(q * q, bool); S = []
        for v in order:
            if not blocked[v]:
                S.append(int(v)); blocked[v] = True; blocked[nb[v]] = True
        if len(S) > len(best): best = S
    return best


@pytest.mark.parametrize("q,kind", [(11, 'std'), (13, 'std')])
def test_an_independent_set_satisfies_every_constraint(q, kind):
    from threepoint import setup
    S = large_independent_set(q, kind)
    n = q * q
    D = build(q, kind, verbose=False, K0=len(S), tri=True, trilocal=True, edge=True)
    cls, lab, zid, gcol = D['cls'], D['lab'], D['zid'], D['gcol']
    n0, Qf, sig, tau, N = setup(q, kind)
    X, Y = np.divmod(np.arange(n), q)
    idx = lambda x, y: (x % q) * q + (y % q)
    maps = []; R = np.eye(2, dtype=int)
    for _ in range(N):
        maps += [R.copy(), (R @ tau) % q]; R = (R @ sig) % q
    perms = [idx(M[0, 0] * X + M[0, 1] * Y, M[1, 0] * X + M[1, 1] * Y) for M in maps]
    inS = np.zeros(n, bool); inS[S] = True; Sa = np.array(S)
    v = np.zeros(D['nv'])
    for c, col in gcol.items():
        Kc = np.nonzero(cls == c)[0]
        v[col] = sum(inS[idx(X[x] + X[Kc], Y[x] + Y[Kc])].sum() for x in S) / (len(S) * len(Kc))
    reps = {}
    for p in range(n * n):
        if zid[lab[p]] >= 0: reps.setdefault(lab[p], p)
    for L_, p in reps.items():
        a, b = divmod(p, n)
        tot = sum((inS[idx(X[Sa] + X[P[a]], Y[Sa] + Y[P[a]])] & inS[idx(X[Sa] + X[P[b]], Y[Sa] + Y[P[b]])]).sum()
                  for P in perms)
        v[zid[L_]] = tot / (len(S) * len(perms))
    assert abs(1 + sum(D['size_c'][c] * v[col] for c, col in gcol.items()) - len(S)) < 1e-9
    for size, ent in D['blocks']:
        assert np.linalg.eigvalsh(block_values(ent, size, v)).min() > -1e-9
    assert max(sum(c * v[i] for i, c in d.items()) - rhs for d, rhs in D['lin']) < 1e-9
    assert (D['Pm'][:, 0] + sum(D['Pm'][:, c] * v[gcol[c]] for c in D['free'])).min() > -1e-9


def test_the_stored_certificates_match_their_checksums():
    """data/threepoint/SHA256SUMS lists every certificate; the slow checks are in test_threepoint_certificates.py"""
    import hashlib
    folder = os.path.join(HERE, '..', 'data', 'threepoint')
    listed = {}
    with open(os.path.join(folder, 'SHA256SUMS')) as f:
        for line in f:
            digest, name = line.split()
            listed[name] = digest
    assert sorted(listed) == sorted(f for f in os.listdir(folder) if f.endswith('.npz'))
    for name, digest in listed.items():
        with open(os.path.join(folder, name), 'rb') as f:
            assert hashlib.sha256(f.read()).hexdigest() == digest, name
