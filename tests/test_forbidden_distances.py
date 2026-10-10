"""Density certificates for planar sets avoiding two distances (papers/forbidden-distances/).

Each certificate in data/forbidden_distances/certificates/ is a dual solution of the linear program of Oliveira
Filho and Vallentin (with clique inequalities for {1, 2} and {1, sqrt 3}), in exact rationals. Three programs check
it: certify_lp2.py, which produced it (Arb balls, Lipschitz sweep), check2.py (60-digit floating point, separate
code) and check3.py (Arb enclosures over whole intervals, exact clique geometry in Q(sqrt 3)). The fast tests run
check3.py on every certificate, check2.py on the small ones and certify_lp2.py on one input; the slow ones run the
rest.
"""
import json
import os
import subprocess
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, '..', 'data', 'forbidden_distances')

# name -> (k, claimed: bound < 1/k, hence chi(R^2, D) >= k + 1)
CERTS = {
    'cl_s3': 6, 's_vs': 6, 's_10': 6, 's_116': 7, 's_4625': 7, 's_40': 8,   # fast
    'cl_2': 6, 's_100': 9,                                                 # slow
}
SLOW = {'cl_2', 's_100'}


def run(script, *args, timeout=1800):
    return subprocess.run([sys.executable, os.path.join(DIR, script), *args], cwd=DIR, capture_output=True,
                          text=True, timeout=timeout)


def exact_bound(name):
    C = json.load(open(os.path.join(DIR, 'certificates', name + '.cert.json')))
    sys.path.insert(0, DIR)
    from shapes import SHAPES
    z0 = Fraction(C['z0'])
    zs = [Fraction(v) for v in C['z']]
    ys = [(Fraction(c['y']), len(SHAPES[c['shape']])) for c in C['cliques']]
    return (z0 + sum(y for y, _ in ys)) / (z0 + sum(zs) + sum(y * m for y, m in ys))


@pytest.mark.parametrize('name', sorted(n for n in CERTS if n not in SLOW))
def test_certificate_checks(name):
    k = CERTS[name]
    assert exact_bound(name) < Fraction(1, k)
    r = run('check2.py', os.path.join('certificates', name + '.cert.json'), str(k))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f'B < 1/{k}: True' in r.stdout


@pytest.mark.parametrize('name', sorted(CERTS))
def test_certificate_check3(name):
    """The rigorous independent check: whole-interval Arb enclosures and exact clique geometry."""
    k = CERTS[name]
    r = run('check3.py', os.path.join('certificates', name + '.cert.json'), str(k))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f'B < 1/{k}: True' in r.stdout


@pytest.mark.slow
@pytest.mark.parametrize('name', sorted(SLOW))
def test_slow_certificate_checks(name):
    k = CERTS[name]
    assert exact_bound(name) < Fraction(1, k)
    r = run('check2.py', os.path.join('certificates', name + '.cert.json'), str(k))
    assert r.returncode == 0, r.stdout + r.stderr
    assert f'B < 1/{k}: True' in r.stdout


def test_certify_reproduces_a_certificate(tmp_path):
    """certify_lp2.py (Arb) rebuilds the stored certificate for {1, 11/6} byte for byte."""
    out = tmp_path / 's_116.cert.json'
    r = run('certify_lp2.py', os.path.join('lp', 's_116.json'), '7', str(out))
    assert r.returncode == 0, r.stdout + r.stderr
    assert 'B < 1/7: True' in r.stdout
    assert out.read_text() == open(os.path.join(DIR, 'certificates', 's_116.cert.json')).read()


@pytest.mark.slow
@pytest.mark.parametrize('name', sorted(CERTS))
def test_certify_reproduces_every_certificate(tmp_path, name):
    out = tmp_path / (name + '.cert.json')
    r = run('certify_lp2.py', os.path.join('lp', name + '.json'), str(CERTS[name]), str(out))
    assert r.returncode == 0, r.stdout + r.stderr
    assert out.read_text() == open(os.path.join(DIR, 'certificates', name + '.cert.json')).read()
