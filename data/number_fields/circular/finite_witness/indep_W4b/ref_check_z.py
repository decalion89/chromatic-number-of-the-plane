#!/usr/bin/env python3
"""Independent referee checker (tasks 1-5).

Claim under test: a finite unit-distance graph H in the plane, with vertex
coordinates in F = Q(sqrt2, sqrt3), has circular chromatic number 4.

  Task 1  exact field arithmetic, distinct points, all-pairs exact distance-1
          test, comparison with the declared edges, float cross-check,
          fixed vertex = origin, degree range
  Task 2  the colouring is proper on the recomputed edges
  Task 3  every listed cycle is simple, has length >= 4 and divisible by 4,
          and uses only edges (closing step included)
  Task 4  own CNF: proper 4-colouring, fixed vertex coloured 0, no listed
          cycle tight (auxiliary tight-arc variables)
  Task 5  kissat -> drat-trim -L -> cake_lpr, with hashes of the formula

Usage:
  python3 ref_check_z.py [--witness FILE] [--no-sat] [--log FILE]
                         [--json-out FILE] [--cnf FILE] [--encoding arc|direct]

Exit code 0: all requested checks passed.  1: rejected.  2: usage error.
Only relative paths are printed.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.normpath(os.path.join(HERE, '..'))
# repository copy: the tools come from the environment variables KISSAT, DRAT_TRIM, CAKE_LPR, else from the PATH
TOOLS = {
    'kissat': os.environ.get('KISSAT') or shutil.which('kissat') or 'kissat',
    'drat-trim': os.environ.get('DRAT_TRIM') or shutil.which('drat-trim') or 'drat-trim',
    'cake_lpr': os.environ.get('CAKE_LPR') or shutil.which('cake_lpr') or 'cake_lpr',
}
EXPECTED_FIELD = 'Q(sqrt2, sqrt3)'
EXPECTED_BASIS = '(1, sqrt2, sqrt3, sqrt6)'
EXPECTED_KEYS = {'field', 'basis', 'denominator', 'points', 'edges',
                 'colouring', 'cycles', 'fixed_vertex'}


def rel(path):
    """Path relative to this folder, for printing (never absolute)."""
    return os.path.relpath(os.path.abspath(path), HERE)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


class Log:
    def __init__(self, path=None):
        self.f = open(path, 'w') if path else None

    def __call__(self, *args):
        s = ' '.join(str(a) for a in args)
        print(s, flush=True)
        if self.f:
            self.f.write(s + '\n')
            self.f.flush()

    def close(self):
        if self.f:
            self.f.close()


# ---------------------------------------------------------------------------
# Field arithmetic in F = Q(sqrt2, sqrt3), basis (1, sqrt2, sqrt3, sqrt6).
# Basis element number p stands for sqrt2^i * sqrt3^j with (i, j) = BASIS[p].
# The table is derived only from sqrt2^2 = 2, sqrt3^2 = 3 and commutativity.
# ---------------------------------------------------------------------------
BASIS = [(0, 0), (1, 0), (0, 1), (1, 1)]
NAMES = ['1', 'sqrt2', 'sqrt3', 'sqrt6']


def derive_table():
    table = {}
    for p, (i1, j1) in enumerate(BASIS):
        for q, (i2, j2) in enumerate(BASIS):
            i, j, coef = i1 + i2, j1 + j2, 1
            if i == 2:          # sqrt2 * sqrt2 = 2
                coef, i = coef * 2, 0
            if j == 2:          # sqrt3 * sqrt3 = 3
                coef, j = coef * 3, 0
            table[p, q] = (coef, BASIS.index((i, j)))
    return table


TABLE = derive_table()


def fmul(x, y):
    """Generic product of two coefficient 4-tuples, driven by TABLE."""
    z = [0, 0, 0, 0]
    for p in range(4):
        for q in range(4):
            c, r = TABLE[p, q]
            z[r] += c * x[p] * y[q]
    return tuple(z)


def fsq(x0, x1, x2, x3):
    """Fast square (x0 + x1 s2 + x2 s3 + x3 s6)^2 written out by hand."""
    return (x0 * x0 + 2 * x1 * x1 + 3 * x2 * x2 + 6 * x3 * x3,
            2 * x0 * x1 + 6 * x2 * x3,
            2 * x0 * x2 + 4 * x1 * x3,
            2 * x0 * x3 + 2 * x1 * x2)


def sqdist_num(p, q):
    """D^2 * |p - q|^2 as a coefficient 4-tuple (exact integers)."""
    a = fsq(p[0] - q[0], p[1] - q[1], p[2] - q[2], p[3] - q[3])
    b = fsq(p[4] - q[4], p[5] - q[5], p[6] - q[6], p[7] - q[7])
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def check_field_arithmetic(log):
    """Print the derived table and prove that fsq agrees with fmul(x, x)."""
    log('  derived multiplication table (row * column):')
    for p in range(4):
        cells = []
        for q in range(4):
            c, r = TABLE[p, q]
            cells.append(('%d*%s' % (c, NAMES[r])) if c != 1 else NAMES[r])
        log('    %-6s | %s' % (NAMES[p], '  '.join('%-9s' % s for s in cells)))
    # Numerical sanity of the table with the positive real roots.
    val = [1.0, math.sqrt(2), math.sqrt(3), math.sqrt(6)]
    worst = 0.0
    for p in range(4):
        for q in range(4):
            c, r = TABLE[p, q]
            worst = max(worst, abs(val[p] * val[q] - c * val[r]))
    log('  table vs real roots: max |error| = %.2e' % worst)
    # Each coordinate of fsq(x) - fmul(x, x) is a polynomial of degree <= 2 in
    # each variable.  Vanishing on the grid S^4 with |S| = 4 > 2 forces it to
    # be the zero polynomial, so agreement on the grid proves the identity.
    grid = [-1, 0, 1, 2]
    ok = all(fsq(*x) == fmul(x, x) for x in itertools.product(grid, repeat=4))
    log('  fast square == generic product on grid {-1,0,1,2}^4 (proves identity):', ok)
    return ok and worst < 1e-12


# ---------------------------------------------------------------------------
# Loading and validation
# ---------------------------------------------------------------------------
def is_int(x):
    return type(x) is int


def load(path, reasons, log):
    try:
        with gzip.open(path, 'rt') as f:
            d = json.load(f)
    except Exception as e:  # noqa: BLE001
        reasons.append(('format', 'cannot read gzipped JSON: %s' % e))
        return None
    if not isinstance(d, dict):
        reasons.append(('format', 'top level is not an object'))
        return None
    missing = EXPECTED_KEYS - set(d)
    if missing:
        reasons.append(('format', 'missing keys %s' % sorted(missing)))
        return None
    extra = set(d) - EXPECTED_KEYS
    if extra:
        log('  note: extra keys ignored:', sorted(extra))
    if d['field'] != EXPECTED_FIELD:
        reasons.append(('format', 'field is %r' % d['field']))
    if d['basis'] != EXPECTED_BASIS:
        reasons.append(('format', 'basis is %r' % d['basis']))
    D = d['denominator']
    if not is_int(D) or D <= 0:
        reasons.append(('format', 'denominator not a positive integer'))
        return None
    P = d['points']
    if not isinstance(P, list) or not all(
            isinstance(p, list) and len(p) == 8 and all(is_int(c) for c in p) for p in P):
        reasons.append(('format', 'points must be lists of 8 integers'))
        return None
    n = len(P)
    if n == 0:
        reasons.append(('format', 'no points'))
        return None
    fv = d['fixed_vertex']
    if not is_int(fv) or not 0 <= fv < n:
        reasons.append(('format', 'fixed_vertex not a valid index'))
        return None
    return d


# ---------------------------------------------------------------------------
# Task 1: points and edges
# ---------------------------------------------------------------------------
def exact_unit_pairs(P, D):
    """All pairs i < j with |p_i - p_j| = 1, by exact integer arithmetic.
    Returns the edge list and the count of pairs whose rational part matches
    D^2 but whose irrational part does not (near misses in the exact sense)."""
    target = (D * D, 0, 0, 0)
    n = len(P)
    edges = []
    rational_only = 0
    for i in range(n):
        p = P[i]
        for j in range(i + 1, n):
            s = sqdist_num(p, P[j])
            if s == target:
                edges.append((i, j))
            elif s[0] == target[0]:
                rational_only += 1
    return edges, rational_only


def float_pass(P, D, exact_edge_set, log):
    """Independent floating-point pass; returns (agrees, info)."""
    import numpy as np
    A = np.array(P, dtype=np.float64)
    r2, r3, r6 = math.sqrt(2.0), math.sqrt(3.0), math.sqrt(6.0)
    X = (A[:, 0] + A[:, 1] * r2 + A[:, 2] * r3 + A[:, 3] * r6) / D
    Y = (A[:, 4] + A[:, 5] * r2 + A[:, 6] * r3 + A[:, 7] * r6) / D
    n = len(P)
    iu, ju = np.triu_indices(n, 1)
    dx = X[iu] - X[ju]
    dy = Y[iu] - Y[ju]
    dist = np.sqrt(dx * dx + dy * dy)
    margin = np.abs(dist - 1.0)
    float_edge = margin < 1e-9
    is_exact = np.zeros(len(iu), dtype=bool)
    # index of pair (i, j), i < j, in the triu ordering
    for (i, j) in exact_edge_set:
        k = i * n - i * (i + 1) // 2 + (j - i - 1)
        is_exact[k] = True
    disagree = int(np.count_nonzero(float_edge != is_exact))
    max_edge = float(margin[is_exact].max()) if is_exact.any() else float('nan')
    ne = ~is_exact
    k_min = int(np.argmin(np.where(ne, margin, np.inf)))
    info = {
        'pairs': int(len(iu)),
        'float_edges': int(np.count_nonzero(float_edge)),
        'disagreements': disagree,
        'max_edge_margin': max_edge,
        'min_nonedge_margin': float(margin[k_min]),
        'min_nonedge_pair': (int(iu[k_min]), int(ju[k_min])),
        'nonedge_pairs_within_1e-12_of_min': int(np.count_nonzero(
            ne & (margin < float(margin[k_min]) + 1e-12))),
        'nonedge_below_1e-3': int(np.count_nonzero(ne & (margin < 1e-3))),
        'nonedge_below_1e-6': int(np.count_nonzero(ne & (margin < 1e-6))),
        'min_dist': float(dist.min()),
    }
    return disagree == 0, info


def decimal_value(c, D, digits=60):
    """(c0 + c1 s2 + c2 s3 + c3 s6) / D^2 to high precision (cross-check only)."""
    getcontext().prec = digits
    s2, s3, s6 = Decimal(2).sqrt(), Decimal(3).sqrt(), Decimal(6).sqrt()
    return (Decimal(c[0]) + Decimal(c[1]) * s2 + Decimal(c[2]) * s3 +
            Decimal(c[3]) * s6) / Decimal(D * D)


# ---------------------------------------------------------------------------
# Task 4: CNF encoding
# ---------------------------------------------------------------------------
def xvar(v, k):
    """Variable 'vertex v has colour k' (k in 0..3)."""
    return 4 * v + k + 1


def build_cnf(n, edges, cycles, fixed, encoding='arc', include_cycles=True):
    """Return (num_vars, clauses, info).

    x(v,k) <-> vertex v has colour k.
      fixed:   x(f,0)
      ALO:     x(v,0) | x(v,1) | x(v,2) | x(v,3)
      AMO:     -x(v,k) | -x(v,l)                     k < l
      edge:    -x(u,k) | -x(v,k)                     every edge uv, every k
    encoding 'arc' (main): one variable y(u,v) per directed arc on a listed
    cycle, forced true when the arc is tight:
      -x(u,k) | -x(v,k+1 mod 4) | y(u,v)             k = 0..3
      -y(a_0) | ... | -y(a_{m-1})                    every listed cycle
    encoding 'direct' (cross-check, needs m = 0 mod 4): for every cycle and
    every start colour k, not all of x(v_t, k+t mod 4):
      -x(v_0,k) | -x(v_1,k+1) | ... | -x(v_{m-1},k+m-1)
    """
    clauses = [[xvar(fixed, 0)]]
    for v in range(n):
        clauses.append([xvar(v, k) for k in range(4)])
        for k in range(4):
            for l in range(k + 1, 4):
                clauses.append([-xvar(v, k), -xvar(v, l)])
    for (u, v) in edges:
        for k in range(4):
            clauses.append([-xvar(u, k), -xvar(v, k)])
    nvars = 4 * n
    info = {'colour_vars': 4 * n, 'arc_vars': 0, 'base_clauses': len(clauses),
            'arc_clauses': 0, 'cycle_clauses': 0}
    if include_cycles:
        if encoding == 'arc':
            arcvar = {}
            for cyc in cycles:
                m = len(cyc)
                for t in range(m):
                    a = (cyc[t], cyc[(t + 1) % m])
                    if a not in arcvar:
                        nvars += 1
                        arcvar[a] = nvars
            for (u, v), y in arcvar.items():
                for k in range(4):
                    clauses.append([-xvar(u, k), -xvar(v, (k + 1) % 4), y])
            info['arc_vars'] = len(arcvar)
            info['arc_clauses'] = 4 * len(arcvar)
            for cyc in cycles:
                m = len(cyc)
                clauses.append([-arcvar[(cyc[t], cyc[(t + 1) % m])] for t in range(m)])
                info['cycle_clauses'] += 1
        elif encoding == 'direct':
            for cyc in cycles:
                m = len(cyc)
                if m % 4:
                    continue  # such a cycle can never be tight
                for k in range(4):
                    clauses.append([-xvar(cyc[t], (k + t) % 4) for t in range(m)])
                    info['cycle_clauses'] += 1
        else:
            raise ValueError(encoding)
    return nvars, clauses, info


def write_cnf(path, nvars, clauses):
    with open(path, 'w') as f:
        f.write('p cnf %d %d\n' % (nvars, len(clauses)))
        for c in clauses:
            f.write(' '.join(map(str, c)) + ' 0\n')


def read_cnf(path):
    with open(path) as f:
        header = f.readline().split()
        assert header[:2] == ['p', 'cnf'], header
        nv, nc = int(header[2]), int(header[3])
        lits = [int(t) for t in f.read().split()]
    clauses, cur = [], []
    for t in lits:
        if t == 0:
            clauses.append(cur)
            cur = []
        else:
            cur.append(t)
    assert not cur
    return nv, nc, clauses


# ---------------------------------------------------------------------------
# Task 5: solver pipeline
# ---------------------------------------------------------------------------
NAME_LINE = re.compile(r'copyright|\(c\)|@|author', re.I)


def sanitize(text):
    """Drop banner lines with personal names, strip absolute paths."""
    out = []
    for line in text.splitlines():
        if NAME_LINE.search(line):
            continue
        for tool, path in TOOLS.items():
            line = line.replace(path, tool)
        line = line.replace(HERE + os.sep, '').replace(HERE, '.')
        line = line.replace(SCRATCH, '..')
        out.append(line)
    return '\n'.join(out) + '\n'


def run_tool(name, args, cwd, log_path, timeout):
    t0 = time.monotonic()
    try:
        proc = subprocess.run([TOOLS[name]] + args, cwd=cwd, stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT, timeout=timeout)
        code, raw = proc.returncode, proc.stdout
    except subprocess.TimeoutExpired as e:
        code, raw = 'timeout', (e.stdout or b'')
    dt = time.monotonic() - t0
    text = sanitize(raw.decode('utf-8', 'replace'))
    with open(log_path, 'w') as f:
        f.write('$ %s %s\n' % (name, ' '.join(args)))
        f.write(text)
        f.write('[exit code %s, wall time %.2f s]\n' % (code, dt))
    return code, text, dt


def min_rotation(c):
    c = list(c)
    return min(tuple(c[s:] + c[:s]) for s in range(len(c)))


def key_lines(text, prefixes):
    return [l for l in text.splitlines() if any(l.startswith(p) for p in prefixes)]


def sat_pipeline(cnf, log, results_dir, timeout, keep_proofs, tag):
    work = os.path.dirname(os.path.abspath(cnf))
    base = os.path.basename(cnf)
    stem = base[:-4] if base.endswith('.cnf') else base
    drat, lrat = stem + '.drat', stem + '.lrat'
    drat_abs, lrat_abs = os.path.join(work, drat), os.path.join(work, lrat)
    hashes = [('before kissat', sha256_file(cnf))]
    res = {}

    code, out, dt = run_tool('kissat', [base, drat], work,
                             os.path.join(results_dir, 'kissat_%s.log' % tag), timeout)
    hashes.append(('after kissat', sha256_file(cnf)))
    s_lines = key_lines(out, ['s '])
    log('  kissat: exit code %s, %.2f s, lines %s' % (code, dt, s_lines))
    res['kissat'] = {'exit': code, 'time': dt, 's_lines': s_lines}
    ok = code == 20 and s_lines == ['s UNSATISFIABLE']
    if os.path.exists(drat_abs):
        res['drat_bytes'] = os.path.getsize(drat_abs)
        log('  DRAT proof size: %d bytes' % res['drat_bytes'])

    if ok:
        code, out, dt = run_tool('drat-trim', [base, drat, '-L', lrat], work,
                                 os.path.join(results_dir, 'drat_trim_%s.log' % tag), timeout)
        hashes.append(('after drat-trim', sha256_file(cnf)))
        s_lines = key_lines(out, ['s '])
        log('  drat-trim: exit code %s, %.2f s, lines %s' % (code, dt, s_lines))
        res['drat-trim'] = {'exit': code, 'time': dt, 's_lines': s_lines}
        ok = 's VERIFIED' in s_lines and os.path.exists(lrat_abs)
        if os.path.exists(lrat_abs):
            res['lrat_bytes'] = os.path.getsize(lrat_abs)
            log('  LRAT proof size: %d bytes' % res['lrat_bytes'])

    if ok:
        code, out, dt = run_tool('cake_lpr', [base, lrat], work,
                                 os.path.join(results_dir, 'cake_lpr_%s.log' % tag), timeout)
        hashes.append(('after cake_lpr', sha256_file(cnf)))
        lines = [l for l in out.splitlines() if l.strip()]
        log('  cake_lpr: exit code %s, %.2f s, output %s' % (code, dt, lines))
        res['cake_lpr'] = {'exit': code, 'time': dt, 'lines': lines}
        ok = code == 0 and 's VERIFIED UNSAT' in lines

    for label, h in hashes:
        log('  formula sha256 %-17s %s' % (label, h))
    same = len(set(h for _, h in hashes)) == 1
    log('  formula unchanged across runs:', same)
    res['hashes'] = hashes
    res['unchanged'] = same
    if not keep_proofs:
        for path in (drat_abs, lrat_abs):
            if os.path.exists(path):
                os.remove(path)  # exact absolute path, no glob
        log('  proof files deleted')
    return ok and same, res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--witness', default=os.path.join(SCRATCH, 'witness_q2_3.json.gz'))   # repository copy: parent folder
    ap.add_argument('--no-sat', action='store_true', help='skip task 5 (solver runs)')
    ap.add_argument('--log', default=None, help='also write this log file')
    ap.add_argument('--json-out', default=None, help='machine-readable summary')
    ap.add_argument('--cnf', default=os.path.join(HERE, 'work', 'ref_z.cnf'))
    ap.add_argument('--encoding', default='arc', choices=['arc', 'direct'])
    ap.add_argument('--results', default=os.path.join(HERE, 'results'))
    ap.add_argument('--timeout', type=float, default=1800.0)
    ap.add_argument('--keep-proofs', action='store_true')
    ap.add_argument('--no-float', action='store_true', help='skip the float pass')
    args = ap.parse_args()

    log = Log(args.log)
    reasons = []
    summary = {'reasons': reasons}

    def reject(code, msg):
        reasons.append((code, msg))
        log('  REJECT-REASON %s: %s' % (code, msg))

    log('Referee check (tasks 1-5)')
    log('witness:', rel(args.witness))
    log('witness sha256:', sha256_file(args.witness))
    log('python', sys.version.split()[0])

    log('\n== Task 1a: field arithmetic ==')
    if not check_field_arithmetic(log):
        reject('internal', 'field arithmetic self-test failed')

    before = len(reasons)
    d = load(args.witness, reasons, log)
    for code, msg in reasons[before:]:
        log('  REJECT-REASON %s: %s' % (code, msg))
    if d is None:
        return finish(log, summary, args, sat_requested=not args.no_sat)
    P, D, fv = d['points'], d['denominator'], d['fixed_vertex']
    n = len(P)
    log('  n = %d points, denominator D = %d, field %r, basis %r'
        % (n, D, d['field'], d['basis']))
    log('  max |coefficient| = %d' % max(abs(c) for p in P for c in p))
    summary['n'] = n

    log('\n== Task 1b: distinct points ==')
    distinct = len(set(map(tuple, P)))
    log('  distinct 8-tuples: %d of %d' % (distinct, n))
    if distinct != n:
        dup = [k for k, c in Counter(map(tuple, P)).items() if c > 1]
        reject('points_not_distinct', '%d repeated tuples, e.g. %s' % (len(dup), list(dup[0])))

    log('\n== Task 1c: exact distance test on all pairs ==')
    t0 = time.monotonic()
    exact_edges, rational_only = exact_unit_pairs(P, D)
    log('  pairs tested: %d (= n(n-1)/2 = %d), time %.1f s'
        % (n * (n - 1) // 2, n * (n - 1) // 2, time.monotonic() - t0))
    log('  unit-distance pairs found: %d' % len(exact_edges))
    log('  pairs with rational part = D^2 but nonzero irrational part: %d' % rational_only)
    exact_set = set(exact_edges)
    summary['exact_edges'] = len(exact_edges)

    log('\n== Task 1d: compare with the declared edges ==')
    E = d['edges']
    declared = set()
    bad = 0
    if not isinstance(E, list):
        reject('edges_malformed', 'edges is not a list')
        E = []
    for e in E:
        if not (isinstance(e, list) and len(e) == 2 and all(is_int(x) for x in e)
                and all(0 <= x < n for x in e) and e[0] != e[1]):
            bad += 1
            continue
        key = (min(e), max(e))
        if key in declared:
            bad += 1
        declared.add(key)
    log('  declared edges: %d entries, %d distinct valid pairs' % (len(E), len(declared)))
    if bad:
        reject('edges_malformed', '%d malformed or repeated edge entries' % bad)
    extra = sorted(declared - exact_set)
    missing = sorted(exact_set - declared)
    log('  declared but not at distance 1: %d' % len(extra))
    log('  at distance 1 but not declared: %d' % len(missing))
    summary['edge_extra'] = extra[:50]
    summary['edge_missing'] = missing[:50]
    if extra or missing:
        reject('edge_mismatch', '%d declared non-unit pairs %s, %d undeclared unit pairs %s'
               % (len(extra), extra[:5], len(missing), missing[:5]))
    else:
        log('  declared edge set == exact unit-distance pairs: True')

    if not args.no_float:
        log('\n== Task 1e: floating-point cross-check ==')
        agree, info = float_pass(P, D, exact_set, log)
        for k, v in info.items():
            log('  %s: %s' % (k, v))
        summary['float'] = info
        i, j = info['min_nonedge_pair']
        c = sqdist_num(P[i], P[j])
        val = decimal_value(c, D)
        log('  closest non-edge pair %s: D^2*|d|^2 = %s + %s*sqrt2 + %s*sqrt3 + %s*sqrt6'
            % ((i, j), c[0], c[1], c[2], c[3]))
        log('  |d|^2 to 60 digits = %s' % val)
        log('  | |d|^2 - 1 | = %.6e' % abs(float(val - 1)))
        if not agree:
            reject('float_crosscheck', 'float pass disagrees with exact pass on %d pairs'
                   % info['disagreements'])

    log('\n== Task 1f: fixed vertex ==')
    log('  fixed vertex %d has coordinates %s' % (fv, P[fv]))
    if any(P[fv]):
        reject('fixed_not_origin', 'fixed vertex %d is not the origin' % fv)
    else:
        log('  fixed vertex is the origin: True')

    adj = [set() for _ in range(n)]
    for (u, v) in exact_edges:
        adj[u].add(v)
        adj[v].add(u)
    deg = [len(a) for a in adj]
    log('\n== Task 1g: degrees (recomputed graph) ==')
    log('  min degree %d, max degree %d, average %.3f, isolated vertices %d'
        % (min(deg), max(deg), sum(deg) / n, deg.count(0)))
    summary['degree_range'] = (min(deg), max(deg))

    log('\n== Task 2: colouring ==')
    col = d['colouring']
    col_ok = (isinstance(col, list) and len(col) == n
              and all(is_int(c) and 0 <= c <= 3 for c in col))
    if not col_ok:
        reject('colouring_malformed', 'colouring must be n integers in {0,1,2,3}')
    else:
        log('  colour counts:', dict(sorted(Counter(col).items())))
        mono = [(u, v) for (u, v) in exact_edges if col[u] == col[v]]
        log('  monochromatic edges: %d' % len(mono))
        if mono:
            reject('colouring_improper', '%d monochromatic edges, e.g. %s' % (len(mono), mono[:5]))
        else:
            log('  colouring is proper on the recomputed edges: True')
        summary['monochromatic'] = mono[:50]

    log('\n== Task 3: cycles ==')
    C = d['cycles']
    cyc_ok = isinstance(C, list) and all(
        isinstance(c, list) and all(is_int(v) and 0 <= v < n for v in c) for c in C)
    if not cyc_ok:
        reject('cycle_malformed', 'cycles must be lists of valid vertex indices')
        C = []
    lengths = Counter(len(c) for c in C)
    log('  cycles listed: %d' % len(C))
    log('  length distribution:', dict(sorted(lengths.items())))
    bad_len = [k for k, c in enumerate(C) if len(c) < 4 or len(c) % 4]
    not_simple = [k for k, c in enumerate(C) if len(set(c)) != len(c)]
    nonedge = []
    for k, c in enumerate(C):
        m = len(c)
        for t in range(m):
            u, v = c[t], c[(t + 1) % m]
            if (min(u, v), max(u, v)) not in exact_set:
                nonedge.append((k, t, m))
                break
    if bad_len:
        reject('cycle_length', '%d cycles with length < 4 or not divisible by 4, e.g. #%d'
               % (len(bad_len), bad_len[0]))
    if not_simple:
        reject('cycle_not_simple', '%d cycles repeat a vertex, e.g. #%d'
               % (len(not_simple), not_simple[0]))
    if nonedge:
        k, t, m = nonedge[0]
        reject('cycle_nonedge', '%d cycles use a non-edge, e.g. cycle #%d step %d of %d%s'
               % (len(nonedge), k, t, m, ' (closing step)' if t == m - 1 else ''))
    summary['cycle_nonedge'] = nonedge[:50]
    summary['cycle_length_bad'] = bad_len[:50]
    summary['cycle_not_simple'] = not_simple[:50]
    if not (bad_len or not_simple or nonedge) and C:
        log('  all cycles simple, length >= 4 and = 0 mod 4, all steps are edges: True')
    canon = set(min_rotation(c) for c in C if c)
    log('  distinct directed cycles up to rotation: %d' % len(canon))
    rev = sum(1 for c in canon if min_rotation(list(reversed(c))) in canon)
    log('  directed cycles whose reverse is also listed: %d' % rev)
    covered = set(v for c in C for v in c)
    log('  vertices on at least one cycle: %d of %d' % (len(covered), n))
    if col_ok and C:
        tight = [k for k, c in enumerate(C)
                 if all((col[c[(t + 1) % len(c)]] - col[c[t]]) % 4 == 1 for t in range(len(c)))]
        log('  cycles tight under the given colouring: %d (lengths %s)'
            % (len(tight), dict(Counter(len(C[k]) for k in tight))))
        summary['tight_given'] = len(tight)

    log('\n== Task 4: CNF encoding (%s) ==' % args.encoding)
    if reasons:
        log('  skipped: the witness was already rejected')
        return finish(log, summary, args, sat_requested=not args.no_sat)
    nvars, clauses, info = build_cnf(n, exact_edges, C, fv, args.encoding)
    os.makedirs(os.path.dirname(os.path.abspath(args.cnf)), exist_ok=True)
    write_cnf(args.cnf, nvars, clauses)
    nv2, nc2, back = read_cnf(args.cnf)
    roundtrip = (nv2 == nvars and nc2 == len(clauses) and back == clauses)
    log('  file: %s' % rel(args.cnf))
    for k, v in info.items():
        log('  %s: %d' % (k, v))
    log('  variables: %d, clauses: %d' % (nvars, len(clauses)))
    log('  clause length distribution:', dict(sorted(Counter(map(len, clauses)).items())))
    log('  re-read file equals generated clauses: %s' % roundtrip)
    h = sha256_file(args.cnf)
    log('  sha256: %s' % h)
    summary['cnf'] = {'vars': nvars, 'clauses': len(clauses), 'sha256': h, **info}
    if not roundtrip:
        reject('internal', 'CNF file round trip failed')
        return finish(log, summary, args, sat_requested=not args.no_sat)

    if args.no_sat:
        log('\n== Task 5: skipped (--no-sat) ==')
        return finish(log, summary, args, sat_requested=False)

    log('\n== Task 5: kissat -> drat-trim -L -> cake_lpr ==')
    os.makedirs(args.results, exist_ok=True)
    ok, res = sat_pipeline(args.cnf, log, args.results, args.timeout, args.keep_proofs,
                           args.encoding)
    summary['sat'] = res
    if not ok:
        reject('sat_pipeline', 'UNSAT was not established and verified')
    return finish(log, summary, args, sat_requested=True)


def finish(log, summary, args, sat_requested):
    reasons = summary['reasons']
    if reasons:
        codes = sorted(set(c for c, _ in reasons))
        log('\nRESULT: REJECT  reasons=%s' % ','.join(codes))
        code = 1
    else:
        scope = 'tasks 1-5' if sat_requested else 'tasks 1-4, solver runs skipped'
        log('\nRESULT: PASS (%s)' % scope)
        code = 0
    summary['result'] = 'PASS' if code == 0 else 'REJECT'
    if args.json_out:
        with open(args.json_out, 'w') as f:
            json.dump(summary, f, indent=1, default=str)
    log.close()
    return code


if __name__ == '__main__':
    sys.exit(main())
