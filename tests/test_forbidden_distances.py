"""Density certificates for planar sets avoiding two distances (papers/forbidden-distances/).

Each certificate in data/forbidden_distances/certificates/ is a dual solution of the linear program of Oliveira
Filho and Vallentin (with clique inequalities for {1, 2} and {1, sqrt 3}), in exact rationals. Three programs check
it: certify_lp2.py, which produced it (Arb balls, Lipschitz sweep), check2.py (60-digit floating point, separate
code) and check3.py (Arb enclosures over whole intervals, exact clique geometry in Q(sqrt 3)). The fast tests run
check3.py on every certificate, check2.py on the small ones and certify_lp2.py on one input; the slow ones run the
rest.

The uniform certificates of data/forbidden_distances/uniform/ (Theorem 5 of the paper) cover whole ranges of d:
tail.py and uniform.py (Taylor enclosures on t-boxes) and check_uniform.py (separate code: Lipschitz enclosures,
boxes in t and d) must accept all of them and reject corrupted ones.
"""
import json
import os
import subprocess
import sys
from fractions import Fraction

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, '..', 'data', 'forbidden_distances')

# name -> k, the claim being bound < 1/k, hence chi(R^2, D) >= k + 1
CERTS = {
    'cl_s3': 6, 's_vs': 6, 's_10': 6, 's_116': 7, 's_4625': 7, 's_40': 8,   # two distances
    'cl_2': 6, 's_100': 9, 's_1000': 11,
    'm3_40': 21, 'm3_100': 26, 'm4_10': 25, 'm4_30': 47,                   # three and four distances
}
SLOW = {'cl_2', 's_100', 's_1000', 'm3_40', 'm3_100', 'm4_10', 'm4_30'}    # check2.py takes minutes or more


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


def _broken(tmp_path, name, mutate):
    C = json.load(open(os.path.join(DIR, 'certificates', name + '.cert.json')))
    mutate(C)
    path = tmp_path / (name + '.broken.cert.json')
    path.write_text(json.dumps(C))
    return str(path)


@pytest.mark.parametrize('name,k,mutate', [
    ('s_116', 7, lambda C: C.update(z0=str(Fraction(C['z0']) - Fraction(1, 100)))),   # g < 0 somewhere
    ('cl_s3', 6, lambda C: C['cliques'][0].update(shape='tri222')),                    # side 2 is not in D
    ('cl_2', 6, lambda C: C['cliques'][0].update(y='-1/1000')),                       # negative weight
])
def test_check3_rejects_broken_certificates(tmp_path, name, k, mutate):
    r = run('check3.py', _broken(tmp_path, name, mutate), str(k))
    assert r.returncode != 0, r.stdout


def test_checkers_fail_when_the_bound_is_too_weak():
    """B = 0.14099... for {1, 11/6} is not below 1/8: both checkers must exit non-zero."""
    for script in ('check3.py', 'check2.py'):
        r = run(script, os.path.join('certificates', 's_116.cert.json'), '8')
        assert r.returncode != 0 and 'B < 1/8: False' in r.stdout, script + r.stdout + r.stderr


@pytest.mark.slow
@pytest.mark.parametrize('name', sorted(SLOW))
def test_slow_certificate_checks(name):
    k = CERTS[name]
    assert exact_bound(name) < Fraction(1, k)
    r = run('check2.py', os.path.join('certificates', name + '.cert.json'), str(k), timeout=None)  # hours for m4_30
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
    r = run('certify_lp2.py', os.path.join('lp', name + '.json'), str(CERTS[name]), str(out), timeout=None)
    assert r.returncode == 0, r.stdout + r.stderr
    assert out.read_text() == open(os.path.join(DIR, 'certificates', name + '.cert.json')).read()


# ----------------------------------------------------------------- uniform certificates (Theorem 5 of the paper)
UDIR = os.path.join(DIR, 'uniform')
TAILS = {7: '163/25', 8: '123/10', 9: '243/10', 10: '529/10', 11: '701/5', 12: '601', 13: '54300'}   # k: D
INTERVALS = {'interval_k7_1.75_1.93': (7, '7/4', '193/100'), 'interval_k8_1.83_1.855': (8, '183/100', '371/200'),
             'interval_k8_12.05_12.3': (8, '241/20', '123/10')}


def urun(script, *args, cwd=UDIR, timeout=1800):
    return subprocess.run([sys.executable, os.path.join(UDIR, script), *args], cwd=cwd, capture_output=True,
                          text=True, timeout=timeout)


def ucert(name):
    return os.path.join('certificates', name + '.json')


def uload(name):
    return json.load(open(os.path.join(UDIR, ucert(name))))


def test_uniform_files_match_the_statement():
    """Theorem 5: chi >= k for d >= D_k (one tail file each) and on the intervals of the three interval files."""
    for k, D in TAILS.items():
        C = uload(f'tail_k{k}')
        z = [Fraction(C[x]) for x in ('z0', 'z1', 'z2')]
        assert C['k'] == k and Fraction(C['D']) == Fraction(D)
        assert z[2] >= 0 and z[0] > 0 and z[0] / sum(z) < Fraction(1, k - 1) and Fraction(C['bound']) == z[0] / sum(z)
    for name, (k, a, b) in INTERVALS.items():
        C = uload(name)
        assert C['k'] == k and Fraction(C['target']) == Fraction(1, k - 1)
        assert Fraction(C['a']) == Fraction(a) and Fraction(C['b']) == Fraction(b)
        ivs = sorted(C['intervals'], key=lambda r: Fraction(r['a']))
        assert Fraction(ivs[0]['a']) == Fraction(a) and Fraction(ivs[-1]['b']) == Fraction(b)
        assert all(Fraction(ivs[i]['b']) == Fraction(ivs[i + 1]['a']) for i in range(len(ivs) - 1))
        for r in ivs:
            z = [Fraction(r[x]) for x in ('z0', 'z1', 'z2')]
            assert Fraction(r['a']) < Fraction(r['b']) and z[0] > 0 and sum(z) > 0
            assert z[0] / sum(z) < Fraction(1, k - 1)
    assert Fraction(INTERVALS['interval_k8_12.05_12.3'][2]) == Fraction(TAILS[8])   # 8 colours from 12.05 on
    assert all(Fraction(TAILS[k]) < Fraction(TAILS[k + 1]) for k in range(7, 13))


def test_uniform_primary_programs_accept():
    r = urun('tail.py', 'check', *[ucert(f'tail_k{k}') for k in TAILS])
    assert r.returncode == 0 and r.stdout.count('CERTIFIED') == len(TAILS), r.stdout + r.stderr
    r = urun('uniform.py', 'check', *[ucert(n) for n in INTERVALS])
    assert r.returncode == 0 and r.stdout.count('CERTIFIED') == len(INTERVALS), r.stdout + r.stderr


def test_uniform_second_checker_accepts():
    r = urun('check_uniform.py', *[ucert(f'tail_k{k}') for k in TAILS], *[ucert(n) for n in INTERVALS])
    assert r.returncode == 0 and r.stdout.count('ACCEPTED') == len(TAILS) + len(INTERVALS), r.stdout + r.stderr


def test_uniform_tail_reproduced(tmp_path):
    """tail.py rebuilds the certificate for 13 colours (d >= 54 300) byte for byte."""
    r = urun('tail.py', '13', '54300', cwd=str(tmp_path))
    assert r.returncode == 0 and 'CERTIFIED' in r.stdout, r.stdout + r.stderr
    assert (tmp_path / 'tail_k13.json').read_text() == open(os.path.join(UDIR, ucert('tail_k13'))).read()


def _first(C, **kw):
    C['intervals'] = sorted(C['intervals'], key=lambda r: Fraction(r['a']))
    C['intervals'][0].update(kw)


UBROKEN = [
    ('tail_k13', lambda C: C.update(z0=str(Fraction(C['z0']) - Fraction(1, 10**5)))),   # h < 0 near t = 3.83
    ('tail_k13', lambda C: C.update(D='50000')),                                        # threshold too low
    ('tail_k13', lambda C: C.update(k=14)),                                             # one colour too many
    ('tail_k7', lambda C: C.update(D='6')),
    ('tail_k9', lambda C: C.update(z2=str(-Fraction(C['z2'])))),                        # z2 < 0: lemma not applicable
    ('interval_k7_1.75_1.93', lambda C: _first(C, z0=str(Fraction(C['intervals'][0]['z0']) - Fraction(1, 100)))),
    ('interval_k8_1.83_1.855', lambda C: _first(C, z2=str(Fraction(C['intervals'][0]['z2']) * Fraction(11, 10)))),
    ('interval_k7_1.75_1.93', lambda C: C['intervals'].pop(5)),                         # a gap
    ('interval_k7_1.75_1.93', lambda C: _first(C, a=C['intervals'][0]['b'], b=C['intervals'][0]['a'])),  # reversed
    ('interval_k7_1.75_1.93', lambda C: (C.update(a='17/10'), _first(C, a='17/10'))),  # below the range of the LP
    ('interval_k8_12.05_12.3', lambda C: C.update(k=9, target='1/8')),                 # one colour too many
]


@pytest.mark.parametrize('i', range(len(UBROKEN)))
def test_uniform_programs_reject_broken_certificates(tmp_path, i):
    name, mutate = UBROKEN[i]
    C = uload(name)
    mutate(C)
    path = tmp_path / (name + '.broken.json')
    path.write_text(json.dumps(C))
    primary = ('tail.py',) if name.startswith('tail') else ('uniform.py',)
    for cmd in (primary + ('check',), ('check_uniform.py',)):
        r = urun(cmd[0], *cmd[1:], str(path))
        assert r.returncode != 0, ' '.join(cmd) + ': ' + r.stdout + r.stderr


@pytest.mark.slow
def test_uniform_programs_reject_a_small_negative_dip(tmp_path):
    """g has minimum -1e-5 on the first interval of interval_k7_1.75_1.93: the floating-point scan of check_uniform.py
    need not see it, and the rigorous parts of both programs must still refuse the certificate."""
    import numpy as np
    from scipy.special import j0
    from scipy.optimize import minimize
    C = uload('interval_k7_1.75_1.93')
    C['intervals'] = sorted(C['intervals'], key=lambda r: Fraction(r['a']))
    r0 = C['intervals'][0]
    a, b = float(Fraction(r0['a'])), float(Fraction(r0['b']))
    z = [float(Fraction(r0[x])) for x in ('z0', 'z1', 'z2')]
    g = lambda v: z[0] + z[1] * j0(v[1]) + z[2] * j0(v[0] * v[1])
    best = None
    for d in np.linspace(a, b, 21):
        t = np.linspace(0, 40, 400001)
        vals = z[0] + z[1] * j0(t) + z[2] * j0(d * t)
        i = int(np.argmin(vals))
        res = minimize(g, [d, t[i]], bounds=[(a, b), (0, 40)], method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-12})
        if best is None or res.fun < best[0]:
            best = (res.fun, res.x)
    r0['z0'] = str(Fraction(r0['z0']) - Fraction(best[0] + 1e-5).limit_denominator(10**12))
    path = tmp_path / 'dip.json'
    path.write_text(json.dumps(C))
    for cmd in (('uniform.py', 'check'), ('check_uniform.py',)):
        r = urun(cmd[0], *cmd[1:], str(path), timeout=None)
        assert r.returncode != 0, ' '.join(cmd) + ': ' + r.stdout + r.stderr
