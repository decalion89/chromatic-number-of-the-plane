"""Density of sets in R^3 avoiding the distance 1 (data/forbidden_distances/m1r3/, papers/unit-distance-density-3d/):
m_1(R^3) <= 0.15158538 for the 27-point configuration X27, below the 0.1532996 of DeCorte, Oliveira Filho and
Vallentin.

The configuration must start with the cuboctahedron and its centre (points q/sqrt2, q in Z^3 with even coordinate sum
and |q|^2 <= 2) and grow by points at distance 1 from at least three earlier points; check_m1r3.py (exact rationals and
Arb) and check_m1r3_indep.py (written separately) must accept the certificate and reject corrupted ones.
"""
import itertools
import json
import os
import shutil
import subprocess
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
MDIR = os.path.join(HERE, '..', 'data', 'forbidden_distances', 'm1r3')
CERT = os.path.join('certificates', 'm1r3_x27.json')
BOUND = '0.15158538'


def mrun(script, *args, timeout=3600):
    return subprocess.run([sys.executable, '-I', os.path.join(MDIR, script), *args], cwd=MDIR, capture_output=True,
                          text=True, timeout=timeout)


def points():
    d = json.load(open(os.path.join(MDIR, 'configurations', 'x27.json')))
    assert d['scale'] == '1/sqrt2'
    return [tuple(Fraction(c) for c in p) for p in d['points']]


def d2(p, q):
    return sum((a - b) ** 2 for a, b in zip(p, q)) / 2


def test_configuration():
    P = points()
    assert len(P) == 27 and len(set(P)) == 27
    fcc = [tuple(Fraction(c) for c in q) for q in itertools.product(range(-1, 2), repeat=3)
           if sum(q) % 2 == 0 and sum(c * c for c in q) <= 2]
    assert P[:13] == fcc
    for k in range(13, 27):
        assert sum(1 for j in range(k) if d2(P[k], P[j]) == 1) >= 3
    edges = [(i, j) for i in range(27) for j in range(i + 1, 27) if d2(P[i], P[j]) == 1]
    assert len(edges) == 82
    sys.path.insert(0, MDIR)
    import conf3
    assert len(conf3.independent_sets(27, edges)) == 44221


def test_certificate_matches_the_configuration_and_beats_the_published_bound():
    C = json.load(open(os.path.join(MDIR, CERT)))
    assert C['side'] == 'M1R3' and C['scale'] == '1/sqrt2'
    assert C['points'] == json.load(open(os.path.join(MDIR, 'configurations', 'x27.json')))['points']
    m, psiT = Fraction(C['m']), Fraction(C['psiT'])
    c = Fraction(C['psi0']) - sum(Fraction(x) for x in C['lam']) - sum(Fraction(e[2]) for e in C['u']) - m
    b = Fraction(BOUND)
    assert c > 0 and c * b * b + m * b - psiT >= 0 and b >= -m / (2 * c)         # the positive root is at most b
    assert b < Fraction(1532996, 10 ** 7)


def test_check_m1r3_accepts():
    r = mrun('check_m1r3.py', CERT)
    assert r.returncode == 0, r.stdout + r.stderr
    assert 'CERTIFIED: every measurable set in R^3 avoiding the distance 1 has upper density at most %s' % BOUND in r.stdout


def test_independent_checker_accepts():
    r = mrun('check_m1r3_indep.py', CERT)
    assert r.returncode == 0 and 'ACCEPTED: upper density <= B <= %s' % BOUND in r.stdout, r.stdout + r.stderr


def _break_congruence(C):
    P = [tuple(Fraction(x) for x in p) for p in C['points']]
    for k, (I, J, v) in enumerate(C['wc']):
        J2 = J[1:] + J[:1]
        if any(d2(P[J2[a] - 1], P[J2[b] - 1]) != d2(P[I[a] - 1], P[I[b] - 1])
               for a in range(len(I)) for b in range(a + 1, len(I))):
            C['wc'][k] = [I, J2, v]
            return
    raise AssertionError('no relation to break')


def _u_on_a_unit_pair(C):
    P = [tuple(Fraction(x) for x in p) for p in C['points']]
    i, j = next((i, j) for i in range(27) for j in range(i + 1, 27) if d2(P[i], P[j]) == 1)
    C['u'][0] = [i + 1, j + 1, C['u'][0][2]]


def _move_a_point(C):
    C['points'][20][0] = str(Fraction(C['points'][20][0]) + Fraction(1, 1000))


def _shift(key, delta):
    def f(C):
        C[key] = str(Fraction(C[key]) + Fraction(delta))
    return f


def _raise_lam(C):
    C['lam'] = [str(Fraction(x) + Fraction(1, 10 ** 5)) for x in C['lam']]


def _index_out_of_range(C):
    C['wc'].append([[1], [28], '1'])


def _psi0_kills_c(C):
    c = Fraction(C['psi0']) - sum(Fraction(x) for x in C['lam']) - sum(Fraction(e[2]) for e in C['u']) - Fraction(C['m'])
    C['psi0'] = str(Fraction(C['psi0']) - c)


CORRUPTIONS = [_break_congruence, _u_on_a_unit_pair, _move_a_point, _shift('psiT', '-1/1000000'),
               _shift('m', '2/100000'), _shift('psi0', '1/100'), _raise_lam, _index_out_of_range, _psi0_kills_c]


@pytest.mark.parametrize('checker', ['check_m1r3.py', 'check_m1r3_indep.py'])
@pytest.mark.parametrize('i', range(len(CORRUPTIONS)))
def test_checkers_reject(tmp_path, i, checker):
    C = json.load(open(os.path.join(MDIR, CERT)))
    CORRUPTIONS[i](C)
    p = tmp_path / 'bad.json'
    json.dump(C, open(p, 'w'))
    r = mrun(checker, str(p))
    assert r.returncode != 0 and 'CERTIFIED' not in r.stdout and 'ACCEPTED' not in r.stdout, r.stdout + r.stderr


@pytest.mark.slow
def test_generator_reproduces_the_certificate(tmp_path):
    """m1cert3.py, from the stored solution of its first program, writes the stored certificate byte for byte"""
    shutil.copy(os.path.join(MDIR, CERT + '.lp.npz'), tmp_path / 'c.json.lp.npz')
    r = subprocess.run([sys.executable, os.path.join(MDIR, 'm1cert3.py'), os.path.join('configurations', 'x27.json'),
                        str(tmp_path / 'c.json'), '--reuse'], cwd=MDIR, capture_output=True, text=True, timeout=3600)
    assert r.returncode == 0, r.stdout + r.stderr
    assert (tmp_path / 'c.json').read_bytes() == open(os.path.join(MDIR, CERT), 'rb').read()
