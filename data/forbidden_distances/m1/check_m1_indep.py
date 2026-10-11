#!/usr/bin/env python3
"""
Independent checker for certificates of the one-distance bound  m_1(R^2) <= psiT  obtained from the
inclusion-exclusion linear program of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki, for a finite
point set X given inside the certificate.

Usage:  python3 check_m1_indep.py CERT.json [--lipschitz] [--series-only] [--workers N] [--debug]
  Last line: "ACCEPTED: m_1(R^2) <= <psiT rounded up to 8 decimals> (<n> points)" or "REJECTED";
  exit code 0 / 1.  Options (off by default): --lipschitz repeats the sweep with a first-derivative
  bound; --series-only uses the power series of J0 for all arguments; --workers N uses N processes;
  --debug shows a Python traceback on internal errors.

Certificate (JSON, numbers are exact rationals written as strings, indices are 1-based integers):
  "side": "M1";  "points": [[[a,b,c,d], [a',b',c',d']], ...]  (x = (a + b sqrt3 + c sqrt11 + d sqrt33)
  + i (a' + b' sqrt3 + c' sqrt11 + d' sqrt33));  "psi0", "psiT";  "lam": one value per point;
  "u": [[i, j, v], ...];  "wc": [[I, J, v], ...].

What is checked (all exactly, in Q(sqrt3, sqrt11), except (P), which is proved rigorously):
  (G) no two points coincide; the unit-distance graph consists of the pairs whose squared distance is
      exactly 1; one lam per point; every u entry has 1 <= i < j <= n, no pair is repeated and no pair
      is at distance 1; every wc entry has |I| = |J| >= 1, all indices in 1..n (also for one-point
      lists), no repeated index within I or within J, and |x_I[a] - x_I[b]|^2 = |x_J[a] - x_J[b]|^2
      for all a < b.
  (A) for EVERY independent set S of the unit-distance graph, including the empty set,
        Q(S) = psiT + sum_{i in S} lam_i + sum_{u: i, j in S} v + sum_{wc} v ([I in S] - [J in S]) >= 0.
      All independent sets are enumerated (their number is cross-checked by the recursion
      I(G) = I(G - v) + I(G - N[v])), and Q is computed for all of them at once by a subset-sum
      transform on the family of independent sets, in exact integers over the common denominator
      (int64 only when the sum of all |terms| is below 2^62, Python integers otherwise); a sample
      is recomputed directly from the definition with exact fractions.
  (P) Phi(t) = psi0 J0(t) - sum_i lam_i - sum_u v J0(|x_i - x_j| t) >= 1 for all t >= 0.
      Integer interval arithmetic (200 fractional bits, outward rounding); J0 enclosed by its power
      series for x <= 30 and by its large-argument expansion (remainder at most sqrt(2/(pi x)) times
      the first omitted term) for x > 30; pi from Machin's formula; distances enclosed from their exact
      squares.  On [0, T]: adaptive steps [a, b] with Phi >= min(Phi(a), Phi(b)) - M2 (b-a)^2/8,
      M2 = |psi0| + sum |v| |x_i - x_j|^2.  For t >= T: Phi >= -sum lam - C/sqrt(t).

Classical facts assumed:
  |J0''| <= 1 and |J0'| <= 1 for real x (Bessel's integral J0(x) = (1/pi) int_0^pi cos(x sin s) ds);
  |J1| <= 1/sqrt(2) (only with --lipschitz);  |J0(x)| <= sqrt(2/(pi x)) for x > 0;
  Hankel's integral for H0^(1), which gives the remainder bound of the large-argument expansion.
Requirements: Python 3.8+ and numpy.
"""
import sys, os, time
from fractions import Fraction as F

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m1_indep_core as core

def round_up_8(q):
    c = -((-q.numerator * 10 ** 8) // q.denominator)          # ceil(q * 10^8)
    sign = '-' if c < 0 else ''
    a, b = divmod(abs(c), 10 ** 8)
    return '%s%d.%08d' % (sign, a, b)

def main():
    args = sys.argv[1:]
    if not args or args[0].startswith('-'):
        print(__doc__); sys.exit(2)
    t0 = time.time()
    path = args[0]
    workers = int(args[args.index('--workers') + 1]) if '--workers' in args else 1
    opts = dict(lipschitz='--lipschitz' in args, workers=workers)
    core.OPTS['series_only'] = '--series-only' in args
    results = {}
    def verdict(name, ok, msg=''):
        results[name] = bool(ok)
        core.log('[%s] %s %s' % ('PASS' if ok else 'FAIL', name, msg))
    def finish(line):
        core.log('\n' + '=' * 78)
        bad = [k for k, v in results.items() if not v]
        core.log('%d checks passed, %d failed%s' % (len(results) - len(bad), len(bad),
                                                   (': ' + '; '.join(bad)) if bad else ''))
        core.log('time: %.1f s' % (time.time() - t0))
        core.log(line)
        sys.exit(0 if line.startswith('ACCEPTED') else 1)

    core.log('certificate: %s' % os.path.basename(path))
    try:
        cert = core.load_cert(path)
    except Exception as ex:
        core.log('certificate unreadable or malformed: %s: %s' % (type(ex).__name__, getattr(ex, 'strerror', None) or ex))
        finish('REJECTED')
    geo = core.Geometry(cert['points'])
    if geo.n == 23:
        core.log('the points are %sthe published 23-point set in its published order'
                 % ('' if cert['points'] == core.X23 else 'NOT '))
    core.selftest_j0()
    t1 = time.time()
    if not core.check_structure(cert, geo, verdict):
        finish('REJECTED')
    core.log('(G) done in %.1f s' % (time.time() - t1))
    t1 = time.time()
    r = core.check_independent_sets(cert, geo)
    if 'nsets' not in r:
        verdict('(A) independent sets', False, r.get('msg', ''))
        finish('REJECTED')
    verdict('(A) %d independent sets enumerated (count cross-checked)' % r['nsets'], r['enum_ok'])
    verdict('(A) Q(S) >= 0 for every independent set (exact)', r['nneg'] == 0 and r['bad'] == 0,
            '(#negative %d, recomputation mismatches %d)' % (r['nneg'], r['bad']))
    core.log('(A) min Q = %s (~%.6e) at S = %s ; Q(empty) = psiT = %s ; %.1f s'
             % (r['minQ'], float(r['minQ']), r['argmin'], r['Qempty'], time.time() - t1))
    t1 = time.time()
    p = core.check_phi(cert, geo, opts, verdict)
    core.log('(P) done in %.1f s' % (time.time() - t1))
    core.log('\nsummary: n = %d points, %d unit-distance pairs, %d independent sets; min Q = %.6e; '
             'smallest certified margin of Phi - 1 = %s; psiT = %s (~%.10f)'
             % (geo.n, len(geo.edges), r['nsets'], float(r['minQ']),
                '%.4e' % float(p['margin']) if p else 'n/a', cert['psiT'], float(cert['psiT'])))
    if p is None or not all(results.values()):
        finish('REJECTED')
    finish('ACCEPTED: m_1(R^2) <= %s (%d points)' % (round_up_8(cert['psiT']), geo.n))

if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as ex:
        if '--debug' in sys.argv:
            raise
        print('error: %s: %s' % (type(ex).__name__, ex))
        print('REJECTED')
        sys.exit(1)
