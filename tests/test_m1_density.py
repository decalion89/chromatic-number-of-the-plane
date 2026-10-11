"""Density of planar sets avoiding the distance 1 (data/forbidden_distances/m1/, papers/unit-distance-density/):
m_1(R^2) <= psiT for the configurations X23 (Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki), X26, X29, X30 and X31.

The configurations must extend each other, with every new point in Q(sqrt3, sqrt11) at distance 1 from at least two
earlier points; check_m1.py (exact arithmetic and Arb) and check_m1_indep.py (written separately: its own field
arithmetic, interval arithmetic and J0) must accept the certificates and reject corrupted ones.
"""
import json
import os
import subprocess
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
MDIR = os.path.join(HERE, '..', 'data', 'forbidden_distances', 'm1')
IEDIR = os.path.join(HERE, '..', 'data', 'forbidden_distances', 'uniform', 'ie')
# name: (points, unit pairs, independent sets, bound on m_1 rounded up to 8 decimals)
CONFIGS = {'x23': (23, 47, 13552, '0.24698967'), 'x26': (26, 60, 28622, '0.24522384'), 'x29': (29, 76, 55610, '0.24414455'),
           'x30': (30, 80, 84893, '0.24379682'), 'x31': (31, 84, 118421, '0.24352462')}
ORDER = ['x23', 'x26', 'x29', 'x30', 'x31']


def mrun(script, *args, timeout=3600):
    return subprocess.run([sys.executable, os.path.join(MDIR, script), *args], cwd=MDIR, capture_output=True, text=True,
                          timeout=timeout)


def cert_path(name):
    return os.path.join('certificates', 'm1_%s.json' % name)


def load_cert(name):
    return json.load(open(os.path.join(MDIR, cert_path(name))))


def exact_points(name):
    sys.path.insert(0, MDIR)
    import exactconf as ec
    pts, _ = ec.load(os.path.join(MDIR, 'configurations', name + '.json'))
    return pts


def test_x23_is_the_configuration_of_ambrus_et_al():
    sys.path.insert(0, IEDIR)
    from kfield import K
    from x23 import P
    assert exact_points('x23') == [(K(re), K(im)) for re, im in P]


@pytest.mark.parametrize('name', sorted(CONFIGS))
def test_configuration(name):
    sys.path.insert(0, MDIR)
    import exactconf as ec
    from kfield import ONE, dist2
    pts = exact_points(name)
    n, ne, na, _ = CONFIGS[name]
    d2, edges = ec.graph(pts)
    assert len(pts) == n and len(edges) == ne and len(ec.independent_sets(n, edges)) == na
    assert len({p for p in pts}) == n
    k = ORDER.index(name); prev = ORDER[k - 1] if k else None
    if prev:
        old = exact_points(prev)
        assert pts[:len(old)] == old
        for k in range(len(old), n):
            assert sum(1 for j in range(k) if dist2(pts[k], pts[j]) == ONE) >= 2


def test_certificates_match_the_configurations_and_improve():
    last = Fraction(1)
    for name in ORDER:
        C = load_cert(name)
        assert C['side'] == 'M1'
        assert C['points'] == json.load(open(os.path.join(MDIR, 'configurations', name + '.json')))['points']
        b = Fraction(C['psiT'])
        assert b < last and b < Fraction(1, 4)
        last = b
    assert last < Fraction(24353, 100000)


@pytest.mark.parametrize('name', sorted(CONFIGS))
def test_check_m1_accepts(name):
    r = mrun('check_m1.py', cert_path(name))
    assert r.returncode == 0 and 'CERTIFIED' in r.stdout, r.stdout + r.stderr
    if CONFIGS[name][3]:
        assert 'at most psiT = %s' % CONFIGS[name][3] in r.stdout


@pytest.mark.parametrize('name', sorted(CONFIGS))
def test_independent_checker_accepts(name):
    r = mrun('check_m1_indep.py', cert_path(name))
    assert r.returncode == 0 and 'ACCEPTED: m_1(R^2) <= ' in r.stdout, r.stdout + r.stderr
    if CONFIGS[name][3]:
        assert 'ACCEPTED: m_1(R^2) <= %s (%d points)' % (CONFIGS[name][3], CONFIGS[name][0]) in r.stdout


def _break_congruence(C):
    """rotate J in the first relation where that breaks the congruence"""
    sys.path.insert(0, MDIR)
    from kfield import K, dist2
    pts = [(K(tuple(Fraction(x) for x in re)), K(tuple(Fraction(x) for x in im))) for re, im in C['points']]
    d = lambda a, b: dist2(pts[a - 1], pts[b - 1])
    for c, (I, J, v) in enumerate(C['wc']):
        J2 = J[1:] + J[:1]
        if any(d(I[a], I[b]) != d(J2[a], J2[b]) for a in range(len(I)) for b in range(a + 1, len(I))):
            C['wc'][c] = [I, J2, v]
            return
    raise AssertionError('no relation to break')


def _one_point_relation_to_24(C):
    """[[1], [24]] with a large weight: in a checker that ignored the range of the indices, [24 in S] would be 0
    and the weight would be free to raise Q"""
    C['wc'].append([[1], [24], '1'])


def _u_on_a_unit_pair(C):
    C['u'][0] = [1, 2, C['u'][0][2]]                 # points 1 and 2 are 0 and 1


def _move_a_point(C):
    C['points'][22][0][0] = str(Fraction(C['points'][22][0][0]) + Fraction(1, 1000))


def _shift(key, delta):
    def m(C):
        C[key] = str(Fraction(C[key]) + Fraction(delta))
    return m


def _raise_lam(C):
    C['lam'] = [str(Fraction(x) + Fraction(1, 10 ** 5)) for x in C['lam']]


def _u_index_0(C):
    C['u'][0][0] = 0


CORRUPTIONS = [_break_congruence, _one_point_relation_to_24, _u_on_a_unit_pair, _move_a_point,
               _shift('psiT', '-1/1000000'), _shift('psi0', '-1/100'), _raise_lam, _u_index_0]


@pytest.mark.parametrize('checker', ['check_m1.py', 'check_m1_indep.py'])
@pytest.mark.parametrize('i', range(len(CORRUPTIONS)))
def test_checkers_reject(tmp_path, i, checker):
    C = load_cert('x23')
    CORRUPTIONS[i](C)
    p = tmp_path / 'bad.json'
    json.dump(C, open(p, 'w'))
    r = mrun(checker, str(p))
    assert r.returncode != 0 and 'CERTIFIED' not in r.stdout and 'ACCEPTED' not in r.stdout, r.stdout + r.stderr


def test_generator_reproduces_the_x23_certificate(tmp_path):
    """m1cert2.py, from the stored solution of its first program, writes the stored certificate byte for byte"""
    import shutil
    shutil.copy(os.path.join(MDIR, 'certificates', 'm1_x23.json.lp.npz'), tmp_path / 'c.json.lp.npz')
    r = mrun('m1cert2.py', os.path.join('configurations', 'x23.json'), str(tmp_path / 'c.json'), '--reuse')
    assert r.returncode == 0, r.stdout + r.stderr
    assert (tmp_path / 'c.json').read_bytes() == open(os.path.join(MDIR, cert_path('x23')), 'rb').read()
