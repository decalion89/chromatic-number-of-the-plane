"""Inclusion-exclusion certificates for two distances (data/forbidden_distances/uniform/ie/, Propositions 19 and 21 of
papers/forbidden-distances/): chi(R^2, {1, d}) >= 15, 16, 17 for every real algebraic d beyond D_15, D_16, D_17.

check_ie.py (exact arithmetic in Q(sqrt3, sqrt11) and Arb) and check_ie_indep.py (written separately: its own interval
arithmetic and J0 series) must accept the certificates and reject corrupted ones; congruences.py must reproduce the
congruence list; the configuration must be the one of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki.
"""
import json
import os
import subprocess
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
IDIR = os.path.join(HERE, '..', 'data', 'forbidden_distances', 'uniform', 'ie')
CERTS = {15: ('ie_k15', '500000'), 16: ('ie_k16', '70000000'), 17: ('ie_k17', '800000000')}     # k: (file, D)


def irun(script, *args, cwd=IDIR, timeout=3600):
    return subprocess.run([sys.executable, os.path.join(IDIR, script), *args], cwd=cwd, capture_output=True, text=True,
                          timeout=timeout)


def path(name):
    return os.path.join('certificates', name + '.json')


def load(name):
    return json.load(open(os.path.join(IDIR, path(name))))


def bound(C, B=None):
    b = Fraction(C['psiT']) + Fraction(C['tau']) * Fraction(C['d0']) ** 2
    if Fraction(C.get('zeta', '0')):
        b -= Fraction(C['zeta']) * (1 - Fraction(B['psiT']))
    return b


def test_configuration_is_x23():
    sys.path.insert(0, IDIR)
    from kfield import ONE, dist2, K
    from x23 import P
    pts = [(K(re), K(im)) for re, im in P]
    d2 = [dist2(pts[i], pts[j]) for i in range(23) for j in range(i + 1, 23)]
    assert len(pts) == 23 and sum(1 for q in d2 if q == ONE) == 47 and len(set(d2)) == 27


def test_congruences_reproduced(tmp_path):
    import shutil
    for f in ('congruences.py', 'x23.py', 'kfield.py'):
        shutil.copy(os.path.join(IDIR, f), tmp_path)
    r = subprocess.run([sys.executable, str(tmp_path / 'congruences.py')], cwd=tmp_path, capture_output=True, text=True,
                       timeout=600)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / 'cong23.json').read_text() == open(os.path.join(IDIR, 'cong23.json')).read()


def test_files_match_the_statement():
    B = load('ie_B')
    for k, (name, D) in CERTS.items():
        C = load(name)
        assert C['side'] == 'A' and Fraction(C['D']) == Fraction(D)
        assert Fraction(C['eta']) >= 0 and Fraction(C['tau']) >= 0 and Fraction(C.get('zeta', '0')) >= 0
        assert bound(C, B) < Fraction(1, k - 1)
        assert (k == 15) == (Fraction(C.get('zeta', '0')) == 0) and (k == 15 or C['B'] == 'ie_B.json')
    assert Fraction(CERTS[15][1]) < Fraction(CERTS[16][1]) < Fraction(CERTS[17][1])


@pytest.mark.parametrize('k', sorted(CERTS))
def test_check_ie_accepts(k):
    r = irun('check_ie.py', path(CERTS[k][0]), str(k))
    assert r.returncode == 0 and 'CERTIFIED' in r.stdout, r.stdout + r.stderr


@pytest.mark.parametrize('k', sorted(CERTS))
def test_independent_checker_accepts(k):
    r = irun('check_ie_indep.py', path(CERTS[k][0]), str(k))
    assert r.returncode == 0 and 'ACCEPTED: chi(R^2, {1, d}) >= %d for every real algebraic d >= %s' % (k, CERTS[k][1]) \
        in r.stdout, r.stdout + r.stderr


def _scale(key, f):
    def m(C):
        C[key] = str(Fraction(C[key]) * Fraction(f))
    return m


def _break_congruence(C):
    """rotate J in the first congruence where that breaks it (a reversed pair, or an equilateral triple, stays congruent)"""
    sys.path.insert(0, IDIR)
    from kfield import K, dist2
    from x23 import P
    pts = [(K(re), K(im)) for re, im in P]
    d = lambda a, b: dist2(pts[a - 1], pts[b - 1])
    for c, (I, J, v) in enumerate(C['wc']):
        J2 = J[1:] + J[:1]
        if any(d(I[a], I[b]) != d(J2[a], J2[b]) for a in range(len(I)) for b in range(a + 1, len(I))):
            C['wc'][c] = [I, J2, v]
            return
    raise AssertionError('no congruence to break')


def _lower_lam(C):
    """every lambda_i lowered by 1e-6: Q drops by at least 1e-6 on every nonempty set, and min Q is about 1e-9"""
    C['lam'] = [str(Fraction(x) - Fraction(1, 10 ** 6)) for x in C['lam']]


def _index_24(C):
    """a one-point 'congruence' [[1], [24]] (X23 has 23 points), paid for by lam_1 and tau: Q and Phi are unchanged
    and the claimed bound falls to 1/140, so a checker that does not range-check the indices would certify 141"""
    d0, tau = Fraction(C['d0']), Fraction(C['tau'])
    w = 2 * d0 * tau * Fraction(9, 10)
    C['lam'][0] = str(Fraction(C['lam'][0]) - w)
    C['wc'].append([[1], [24], str(w)])
    C['tau'] = str(tau - w / (2 * d0))


def _u_index_0(C):
    C['u'][0][0] = 0


CORRUPTIONS = [
    ('ie_k15', _index_24, None),
    ('ie_k15', _u_index_0, None),
    ('ie_k15', _scale('eta', '6/5'), None),
    ('ie_k15', _lower_lam, None),
    ('ie_k15', _break_congruence, None),
    ('ie_k15', lambda C: C.update(D=str(Fraction(C['D']) / 10)), None),
    ('ie_k17', _scale('zeta', '51/50'), None),
    ('ie_k17', lambda C: None, lambda B: B.update(pi=[[i, j, str(Fraction(v) * Fraction(9, 10))] for i, j, v in B['pi']])),
    ('ie_k17', lambda C: None, lambda B: B['pi'][0].__setitem__(1, 24)),
]


@pytest.mark.parametrize('checker', ['check_ie.py', 'check_ie_indep.py'])
@pytest.mark.parametrize('i', range(len(CORRUPTIONS)))
def test_checkers_reject(tmp_path, i, checker):
    name, m, mb = CORRUPTIONS[i]
    C = load(name); B = load('ie_B')
    m(C)
    if mb:
        mb(B)
    d = tmp_path / 'c'; d.mkdir()
    json.dump(B, open(d / 'ie_B.json', 'w')); json.dump(C, open(d / (name + '.json'), 'w'))
    k = {'ie_k15': 15, 'ie_k16': 16, 'ie_k17': 17}[name]
    r = irun(checker, str(d / (name + '.json')), str(k))
    assert r.returncode != 0 and 'CERTIFIED' not in r.stdout and 'ACCEPTED' not in r.stdout, r.stdout + r.stderr
