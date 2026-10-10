#!/usr/bin/env python3
"""
Independent checker for certificates of  chi(R^2, {1, d}) >= k  for all real algebraic d >= D,
built on the 23-point set of Ambrus, Csiszarik, Matolcsi, Varga and Zsamboki.

Usage:  python3 check_ie_indep.py CERT.json [k] [--lipschitz] [--ab-crosscheck] [--series-only] [--workers N]

  With k: the bound must satisfy  bound < 1/(k-1).  Without k: the largest such k is reported.
  The last line says ACCEPTED (with the certified statement) or REJECTED; exit code 0 / 1.
  Options (all off by default): --lipschitz repeats every sweep with a first-derivative bound;
  --ab-crosscheck recomputes all 2^23 values of (AB) a second way; --series-only uses the power
  series of J0 for all arguments (no large-argument expansion); --workers N uses N processes;
  --B FILE uses another side-B certificate than the one named in a side-A certificate;
  --debug shows a Python traceback instead of a one-line error message.

Certificates (JSON; numbers are exact rationals written as strings):
  side "B": psi0, psiT, pi, lam (23 values), u = [[i, j, v], ...], wc = [[I, J, v], ...]; pi is one
    rational, or a list [[i, j, pi_e], ...] with exactly one entry for each of the 47 unit pairs (i < j).
  side "A": D, d0, mu, psi0, psiT, eta, tau, lam, u, wc, and optionally zeta with "B" = the file name
    (relative to the A file) of a side-B certificate, which is then checked completely first.
    Without a B file, zeta must be absent or 0.

What is checked (points are 1-based; r_ij = |x_i - x_j|; e(y) = min(mu, sqrt(2/(pi y)))):
  Points: the embedded exact coordinates (in Q(sqrt3, sqrt11, i)) give 47 unit-distance pairs,
    27 distinct distances and the published multiset of all 253 distances.
  Both sides, exactly: 23 lam values; every u pair is a non-unit pair; every wc entry has I, J of the
    same length with distinct indices and |x_I[a] - x_I[b]| = |x_J[a] - x_J[b]| for all a < b.
  Side B:
    (GB) every pi_e >= 0.
    (AB) for EVERY subset T of the 23 points:  psiT + sum_{i in T} lam_i + sum_{u inside T} v
         + sum_{unit pairs e inside T} pi_e + sum_{wc} v ([I in T] - [J in T]) >= 0          (exact).
    (PB) for all t >= 0:  psi0 J0(t) - sum lam - sum_u v J0(r_ij t) >= 1.
    Passed on: psiT_B and L_B = max(0, psi0 + sum_e pi_e).
  Side A:
    (G) eta, tau, zeta >= 0, 0 < d0 < 1, D > 0, mu >= 0.4028.
    (Z) psi0 - sum lam - sum u + eta + tau (2 d0 - 1) >= 1                                    (exact).
    (A) for every nonempty independent set S of the unit-distance graph:
        psiT + sum_{i in S} lam_i + sum_{u inside S} v + sum_{wc} v ([I in S] - [J in S]) >= 0,
        and psiT >= zeta                                                                      (exact).
    (P) for all t > 0:  psi0 J0(t) - sum lam - sum_u v J0(r_ij t) - eta e(Dt) + 2 d0 tau
                        - zeta L_B e(Dt) (529 - H(t)) >= 1,   H(t) = sum_{i,j} J0(r_ij t);
        and J0(x) >= -0.4028 for all x >= 0 (so that -J0(y) <= e(y) for y > 0).
    bound = psiT + tau d0^2 - zeta (1 - psiT_B)   (exact).

How (P) and (PB) are proved: all arithmetic is on integers (fixed point with 200 fractional bits,
every operation rounded outward).  J0 is enclosed by its power series for x <= 30 and by its
large-argument expansion for x > 30, with remainder at most sqrt(2/(pi x)) times the first omitted
term; pi comes from Machin's formula.  On [0, T] the function is checked on an adaptive grid: on each
step [a, b] it is at least the smaller endpoint value minus M2 (b-a)^2/8, M2 bounding the second
derivative, and e(Dt) <= e(Da) there.  For t >= T it is at least its limit minus C/sqrt(t).

Classical facts assumed:
  |J0| <= 1 and |J0''| <= 1 for real x (Bessel's integral J0(x) = (1/pi) int_0^pi cos(x sin s) ds);
  |J1| <= 1/sqrt(2) for real x (from J0^2 + 2 sum_n Jn^2 = 1; used only with --lipschitz);
  |J0(x)| <= sqrt(2/(pi x)) for x > 0;
  Hankel's integral for H0^(1), which gives the remainder bound of the large-argument expansion;
  H(t) >= 0 (H(t) is the average over directions w of |sum_k exp(i t <x_k, w>)|^2).
Requirements: Python 3.8+; numpy for the check over all 2^23 subsets of side-B certificates.
"""
import sys, os, time
from fractions import Fraction as F

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ie_indep_core as core

def main():
    args = sys.argv[1:]
    if not args or args[0].startswith('-'):
        print(__doc__); sys.exit(2)
    t0 = time.time()
    path = args[0]
    def opt(name, conv=str, default=None):
        return conv(args[args.index(name) + 1]) if name in args else default
    flags = {a for a in args if a.startswith('--')}
    vals = [a for i, a in enumerate(args[1:], 1) if not a.startswith('--') and args[i - 1] not in ('--workers', '--B')]
    if vals and not (vals[0].isdigit() and int(vals[0]) >= 2):
        print('k must be an integer >= 2'); sys.exit(2)
    k = int(vals[0]) if vals else None
    opts = dict(lipschitz='--lipschitz' in flags, ab_crosscheck='--ab-crosscheck' in flags,
                workers=opt('--workers', int, 1))
    core.OPTS['series_only'] = '--series-only' in flags
    results = {}
    def verdict(name, ok, msg=''):
        results[name] = bool(ok)
        core.log('[%s] %s %s' % ('PASS' if ok else 'FAIL', name, msg))
    def finish(line):
        core.log('\n' + '=' * 78)
        bad = [n for n, ok in results.items() if not ok]
        core.log('%d checks passed, %d failed%s' % (len(results) - len(bad), len(bad),
                                                   (': ' + '; '.join(bad)) if bad else ''))
        core.log('[%.1f s] %s' % (time.time() - t0, line))
        sys.exit(0 if line.startswith('ACCEPTED') else 1)

    core.log('certificate: %s%s' % (os.path.basename(path), '' if k is None else ', target k = %d' % k))
    try:
        cert = core.load_cert(path)
    except Exception as ex:
        finish('REJECTED: certificate unreadable or malformed (%s: %s)'
               % (type(ex).__name__, getattr(ex, 'strerror', None) or ex))
    geo = core.Geometry()
    core.log('points: 47 unit pairs, 27 distinct distances, distance multiset as published (asserted)')
    core.selftest_j0()

    if cert['kind'] == 'B':
        if k is not None:
            core.log('(k is not used for a side-B certificate)')
        Bres = core.run_B(cert, geo, opts, verdict)
        if Bres is not None and all(results.values()):
            finish('ACCEPTED: side-B certificate valid; psiT_B = %s, L_B = %s' % (Bres['psiT'], Bres['L']))
        finish('REJECTED')

    Bres = None
    if cert['Bfile'] is not None:
        bpath = opt('--B') or os.path.join(os.path.dirname(os.path.abspath(path)), cert['Bfile'])
        core.log('B certificate: %s' % os.path.basename(bpath))
        try:
            certB = core.load_cert(bpath)
            if certB['kind'] != 'B':
                raise ValueError('the referenced certificate is not a side-B certificate')
        except Exception as ex:
            finish('REJECTED: B certificate unreadable or malformed (%s: %s)'
                   % (type(ex).__name__, getattr(ex, 'strerror', None) or ex))
        Bres = core.run_B(certB, geo, opts, verdict)
    res = core.run_A(cert, Bres, geo, opts, verdict)
    if res is None or not all(results.values()):
        finish('REJECTED')
    bound, D = res['bound'], res['D']
    core.log('\nbound = psiT + tau d0^2 - zeta (1 - psiT_B) = %s (~%.12f)' % (bound, float(bound)))
    if bound <= 0:
        finish('REJECTED: the bound is not positive')
    kmax = -((-bound.denominator) // bound.numerator)          # largest k with bound < 1/(k-1) is ceil(1/bound)
    if kmax < 2:
        finish('REJECTED: the bound is >= 1')
    core.log('largest k with bound < 1/(k-1): %d   (1/(k-1) = %.9f; 1/k = %.9f)' % (kmax, 1 / (kmax - 1), 1 / kmax))
    if k is not None:
        if not bound < F(1, k - 1):
            finish('REJECTED: bound %.9f is not < 1/(k-1) = %.9f for k = %d' % (float(bound), 1 / (k - 1), k))
        kk = k
    else:
        kk = kmax
    finish('ACCEPTED: chi(R^2, {1, d}) >= %d for every real algebraic d >= %s' % (kk, D))

if __name__ == '__main__':
    try:
        main()
    except SystemExit:
        raise
    except Exception as ex:                      # report without a traceback (use --debug to see one)
        if '--debug' in sys.argv:
            raise
        print('error: %s: %s' % (type(ex).__name__, ex))
        print('REJECTED')
        sys.exit(1)
