"""An independent check of the three-point certificates (data/threepoint/*.npz).

scripts/threepoint_verify_indep.py does not import scripts/threepoint.py or scripts/threepoint_verify.py. It
rebuilds the programme of notes/local_colourings.md, section 14, from its definitions: the orbits of ordered
pairs under the isometries and the relabellings of a triangle, the matrices M1 and M0, and the triangle and
Delsarte inequalities, over full n x n matrices. From a certificate it takes only the dual multipliers and the
labelling that says what each multiplies. It proves the dual blocks positive definite in exact arithmetic and
bounds alpha by an exact rational. These tests check it on its own terms:
- on small planes, the true z of an actual independent set satisfies every constraint it rebuilds, so the
  constraints are valid and the orbits are not too coarse;
- misreading a certificate's labelling destroys the bound, so the check is not vacuous;
- the stored certificates give the claimed bounds: inert13 and inert29 in CI, all eight marked slow.
"""
import os, sys
from fractions import Fraction
import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path.append(os.path.join(HERE, '..', 'scripts'))   # appended, so that no script can shadow an installed package

pytest.importorskip('flint')
from threepoint_verify_indep import Plane, pair_classes, Check   # noqa: E402

CERTS = os.path.join(HERE, '..', 'data', 'threepoint')

# certificate: (q, the lower bound on chi it proves), as in tests/test_threepoint_certificates.py
CLAIMS = {'inert13.npz': (13, 5), 'inert29.npz': (29, 6), 'std37.npz': (37, 6), 'inert37.npz': (37, 6),
          'std41.npz': (41, 6), 'inert41.npz': (41, 6), 'std43.npz': (43, 6), 'std47.npz': (47, 6)}


def greedy_independent(P, rng):
    """a random maximal independent set of the plane P, the largest of 20 greedy tries"""
    n = P.n
    U = np.nonzero(P.unit)[0]
    nbr = [P.idx(P.X[p] + P.X[U], P.Y[p] + P.Y[U]) for p in range(n)]
    best = None
    for _ in range(20):
        inS = np.zeros(n, bool)
        blocked = np.zeros(n, int)
        for p in rng.permutation(n):
            if blocked[p] == 0:
                inS[p] = True
                blocked[nbr[p]] += 1
        if best is None or inS.sum() > best.sum():
            best = inS.copy()
    S = np.nonzero(best)[0]
    assert not any(np.any(best[nbr[a]]) for a in S)     # independent
    return S


def true_z(P, S):
    """Zint[a, b] = sum over R in O(Q) of #{p in S : p + Ra, p + Rb in S}; z = Zint / (|O(Q)| |S|)"""
    n = P.n
    inS = np.zeros(n, bool)
    inS[S] = True
    Ind = np.array([inS[P.idx(P.X[p] + P.X, P.Y[p] + P.Y)] for p in S], dtype=np.int64)   # [p + w in S]
    K = Ind.T @ Ind
    Zint = np.zeros((n, n), dtype=np.int64)
    for R in P.OQ:
        img = P.image(R)
        Zint += K[np.ix_(img, img)]
    return Zint, len(P.OQ) * len(S)


@pytest.mark.parametrize("q,kind,n0", [(7, 'std', 1), (11, 'inert', 1), (13, 'inert', 2)])
def test_true_z_of_independent_sets_satisfies_the_rebuilt_constraints(q, kind, n0):
    rng = np.random.default_rng(q)
    P = Plane(q, kind, n0)
    n = P.n
    S = greedy_independent(P, rng)
    Zint, den = true_z(P, S)
    z = Zint / den
    uniq, inv = pair_classes(P)
    zi = Zint.ravel()
    mx = np.full(len(uniq), -1, dtype=np.int64)
    mn = np.full(len(uniq), 1 << 62, dtype=np.int64)
    np.maximum.at(mx, inv, zi)
    np.minimum.at(mn, inv, zi)
    assert np.array_equal(mx, mn)                                            # z is constant on every orbit
    ra, rb = uniq // n, uniq % n
    unit_sided = P.unit[ra] | P.unit[rb] | P.unit[P.diff(rb, ra)]
    assert Zint[0, 0] == den and np.all(mx[unit_sided] == 0)
    assert np.all((zi >= 0) & (zi <= den))
    assert Zint[0].sum() == len(S) * den                                     # sum_w z(0, w) = |S|
    D = z[0][P.diff(np.arange(n)[None, :], np.arange(n)[:, None])]            # D[a, b] = z(0, b - a)
    assert np.linalg.eigvalsh(z).min() > -1e-9                               # M1 is PSD
    assert np.linalg.eigvalsh((D - z)[1:, 1:]).min() > -1e-9                 # M0 is PSD
    degen = (ra == 0) | (rb == 0) | (ra == rb)
    for t in np.nonzero(~degen & ~unit_sided)[0]:                            # triangle inequalities, exactly
        a, b = int(ra[t]), int(rb[t])
        g = [Zint[0, a], Zint[0, b], Zint[0, P.diff(b, a)]]
        zt = Zint[a, b]
        assert min(gi - zt for gi in g) >= 0
        assert min(den - g[i] - g[j] + zt for i, j in ((0, 1), (0, 2), (1, 2))) >= 0
    for xi in range(1, n):                                                   # Delsarte
        a_, b_ = divmod(xi, q)
        assert np.sum(z[0] * np.cos(2 * np.pi * ((a_ * P.X + b_ * P.Y) % q) / q)) > -1e-9


def test_misreading_the_labelling_destroys_the_bound():
    ch = Check(os.path.join(CERTS, 'inert13.npz'))
    right = ch.bound_float(ch.psd_coef_float())[0]
    assert 42 < right < 43                                # alpha(G_13) <= 42 < 169/4
    shift = ch.shift
    ch.shift = [-s for s in shift]                        # the half-shifts of the basis with the wrong sign
    assert ch.bound_float(ch.psd_coef_float())[0] > 2 * right
    ch.shift = shift
    cert_blocks = ch.cert_blocks                          # block k of M1 read as block k + 1
    N2 = ch.N // 2 + 1
    perm = list(range(len(ch.blocks)))
    for i in range(1, N2 - 1):
        perm[i] = i + 1
    perm[N2 - 1] = 1
    ch.cert_blocks = [cert_blocks[perm[i]] for i in range(len(ch.blocks))]
    assert ch.bound_float(ch.psd_coef_float())[0] > 2 * right
    ch.cert_blocks = cert_blocks
    assert ch.bound_float(ch.psd_coef_float())[0] == right


def check_certificate(name):
    q, chi = CLAIMS[name]
    ch = Check(os.path.join(CERTS, name))
    assert ch.q == q
    bound = ch.bound_exact(ch.psd_coef_exact(chunk=500))[0]
    assert isinstance(bound, Fraction)
    alpha = bound.numerator // bound.denominator          # alpha <= floor of the exact rational bound
    assert (chi - 1) * alpha < q * q                      # so chi(plane) >= q^2 / alpha > chi - 1


@pytest.mark.parametrize("name", ['inert13.npz', 'inert29.npz'])
def test_small_certificates(name):
    check_certificate(name)


@pytest.mark.slow
@pytest.mark.parametrize("name", sorted(set(CLAIMS) - {'inert13.npz', 'inert29.npz'}))
def test_other_certificates(name):
    check_certificate(name)
