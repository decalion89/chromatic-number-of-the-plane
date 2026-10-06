#!/usr/bin/env python3
"""Sanity tests (task 6) for ref_check_z.py.

Part A  corrupted copies of the witness.  Each copy is checked by
        `ref_check_z.py --no-sat` in a separate process; it must be rejected
        and the set of reason codes must be exactly the expected set (or
        contain it, where the corruption has unavoidable side effects).
        Expectations about geometry are computed here with floats, not with
        the checker's exact code.
Part B  encoding tests with kissat (no proofs needed, these are SAT runs):
        B1  formula without cycle clauses is satisfiable; every model (ten
            runs on randomly renamed copies) decodes to a proper 4-colouring
            with the fixed vertex 0, and has a tight listed cycle.
        B2  the given colouring, rotated so the fixed vertex has colour 0,
            violates exactly the cycle clauses of its tight cycles.
        B3  dropping the cycles that are tight under the given colouring
            makes the formula satisfiable; the model has no tight listed cycle.

Usage: python3 ref_tests_z.py [--log FILE]
Exit code 0 when every test passes.  Only relative paths are printed.
"""
import argparse
import copy
import gzip
import importlib.util
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ref_check_z', os.path.join(HERE, 'ref_check_z.py'))
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

# repository copy: the corrupted files go to a temporary folder (or $REF_TESTDIR); the witness is in the parent folder
TESTDIR = os.environ.get('REF_TESTDIR') or tempfile.mkdtemp(prefix='ref_tests_z_')
WITNESS = os.path.join(os.path.dirname(HERE), 'witness_q2_3.json.gz')
R2, R3, R6 = math.sqrt(2.0), math.sqrt(3.0), math.sqrt(6.0)


def fcoord(p, D):
    return ((p[0] + p[1] * R2 + p[2] * R3 + p[3] * R6) / D,
            (p[4] + p[5] * R2 + p[6] * R3 + p[7] * R6) / D)


def float_unit_neighbours(P, D, v, tol=1e-9):
    x, y = fcoord(P[v], D)
    out = []
    for w, q in enumerate(P):
        if w != v:
            a, b = fcoord(q, D)
            if abs(math.hypot(x - a, y - b) - 1.0) < tol:
                out.append(w)
    return out


def tight(col, c):
    m = len(c)
    return all((col[c[(t + 1) % m]] - col[c[t]]) % 4 == 1 for t in range(m))


class Runner:
    def __init__(self, log):
        self.log = log
        self.failures = []
        self.count = 0

    def check(self, name, cond, detail=''):
        self.count += 1
        self.log('  [%s] %s %s' % ('ok' if cond else 'FAIL', name, detail))
        if not cond:
            self.failures.append(name)


def run_checker(path, json_out, cnf_out):
    t0 = time.monotonic()
    proc = subprocess.run([sys.executable, os.path.join(HERE, 'ref_check_z.py'),
                           '--witness', path, '--no-sat', '--json-out', json_out,
                           '--cnf', cnf_out],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    dt = time.monotonic() - t0
    try:
        with open(json_out) as f:
            summary = json.load(f)
    except OSError:
        summary = {'reasons': [['checker_crashed', proc.stdout.decode()[-500:]]],
                   'result': 'CRASH'}
    return proc.returncode, summary, proc.stdout.decode(), dt


def kissat_solve(cnf, seed=0, timeout=600):
    """Run kissat without a proof; return (status, model dict var->bool)."""
    proc = subprocess.run([rc.TOOLS['kissat'], '--seed=%d' % seed, os.path.basename(cnf)],
                          cwd=os.path.dirname(cnf), stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=timeout)
    out = proc.stdout.decode()
    status = [l for l in out.splitlines() if l.startswith('s ')]
    model = {}
    for l in out.splitlines():
        if l.startswith('v '):
            for t in l[2:].split():
                lit = int(t)
                if lit:
                    model[abs(lit)] = lit > 0
    return proc.returncode, status, model


def decode(model, n):
    """Colour per vertex; None when a vertex does not have exactly one colour."""
    col = []
    for v in range(n):
        ks = [k for k in range(4) if model.get(rc.xvar(v, k), False)]
        if len(ks) != 1:
            return None
        col.append(ks[0])
    return col


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    log = rc.Log(args.log)
    R = Runner(log)
    os.makedirs(TESTDIR, exist_ok=True)

    with gzip.open(WITNESS, 'rt') as f:
        W = json.load(f)
    P, D, fv, col, C = W['points'], W['denominator'], W['fixed_vertex'], W['colouring'], W['cycles']
    n = len(P)
    edges, _ = rc.exact_unit_pairs(P, D)
    eset = set(edges)
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    on_cycle = set(v for c in C for v in c)
    log('Sanity tests (task 6)')
    log('witness: %s, n = %d, exact edges = %d' % (rc.rel(WITNESS), n, len(edges)))

    # ------------------------------------------------------------------
    # Part A: corrupted witnesses
    # ------------------------------------------------------------------
    cases = []

    def case(name, desc, mutate, expected, exact=True, extra_check=None):
        cases.append((name, desc, mutate, expected, exact, extra_check))

    # A0 original passes
    case('original', 'unchanged witness', lambda w: None, set())

    # A1 moved point on a cycle (not the fixed vertex): x += 1/36
    v_on = next(v for v in C[0] if v != fv)

    def move(w, v, da0=1, da1=0, axis=0):
        # axis 0: x-coordinate (a0, a1); axis 4: y-coordinate (b0, b1)
        w['points'][v][axis] += da0
        w['points'][v][axis + 1] += da1

    def extra_pairs_touch(v):
        def chk(s):
            pairs = [tuple(p) for p in s['edge_extra'] + s['edge_missing']]
            return bool(pairs) and all(v in p for p in pairs), 'mismatched pairs all contain %d' % v
        return chk

    case('moved_point_on_cycle',
         'vertex %d (on %d cycles) moved by +1/36 in x' % (v_on, sum(v_on in c for c in C)),
         lambda w: move(w, v_on), {'edge_mismatch', 'cycle_nonedge'},
         extra_check=extra_pairs_touch(v_on))

    # A2 moved point on no cycle (degree >= 1)
    off = sorted(v for v in range(n) if v not in on_cycle)
    v_off = next(v for v in off if adj[v] and v != fv)
    case('moved_point_off_cycle',
         'vertex %d (on no cycle, degree %d) moved by +1/36 in x' % (v_off, len(adj[v_off])),
         lambda w: move(w, v_off), {'edge_mismatch'}, extra_check=extra_pairs_touch(v_off))

    # A3 tiny move by (19601 - 13860 sqrt2)/36, about 7.1e-7, along the axis on
    # which every neighbour has a large offset (so d^2 changes by > 1e-8 and the
    # float view sees it too; the exact test does not need this)
    tiny = (19601 - 13860 * R2) / D
    xo, yo = fcoord(P[v_off], D)
    offs = {0: min(abs(fcoord(P[w], D)[0] - xo) for w in adj[v_off]),
            4: min(abs(fcoord(P[w], D)[1] - yo) for w in adj[v_off])}
    ax = max(offs, key=offs.get)
    case('tiny_move_off_cycle',
         'vertex %d moved by (19601-13860*sqrt2)/36 = %.3e along %s (min neighbour offset %.3f)'
         % (v_off, tiny, 'x' if ax == 0 else 'y', offs[ax]),
         lambda w: move(w, v_off, 19601, -13860, ax), {'edge_mismatch'},
         extra_check=extra_pairs_touch(v_off))

    # A4 improper colouring
    u0, w0 = next((u, w) for (u, w) in edges if u != fv and w != fv)

    def recolour(w):
        w['colouring'][u0] = w['colouring'][w0]
    case('improper_colouring', 'colour of %d set to colour of neighbour %d' % (u0, w0),
         recolour, {'colouring_improper'},
         extra_check=lambda s: ([u0, w0] in s['monochromatic'] or [w0, u0] in s['monochromatic'],
                                'reported monochromatic edge contains %d-%d' % (u0, w0)))

    # A5 cycle through a non-edge (inner step)
    k4 = next(k for k, c in enumerate(C) if len(c) == 4)
    c4 = C[k4]
    w_bad = next(x for x in range(n) if x not in c4 and x not in adj[c4[1]])

    def bad_inner(w):
        w['cycles'][k4] = [c4[0], c4[1], w_bad, c4[3]]
    case('cycle_nonedge_inner', 'cycle #%d: vertex %d replaced by non-neighbour %d of %d'
         % (k4, c4[2], w_bad, c4[1]), bad_inner, {'cycle_nonedge'},
         extra_check=lambda s: (s['cycle_nonedge'][0][:2] == [k4, 1],
                                'first bad step reported as cycle #%d step 1' % k4))

    # A6 cycle whose closing step is the only non-edge
    k8 = c8 = None
    for k, c in enumerate(C):
        if len(c) == 8:
            for r in range(8):
                q = c[r:] + c[:r]
                if (min(q[3], q[0]), max(q[3], q[0])) not in eset:
                    k8, c8 = k, q[:4]
                    break
        if k8 is not None:
            break

    def bad_closing(w):
        w['cycles'][k8] = list(c8)
    case('cycle_nonedge_closing', 'cycle #%d replaced by path %s whose closing step is a non-edge'
         % (k8, c8), bad_closing, {'cycle_nonedge'},
         extra_check=lambda s: (s['cycle_nonedge'][0] == [k8, 3, 4],
                                'reported as cycle #%d step 3 of 4 (closing)' % k8))

    # A7 non-simple closed walk a,b,a,b
    a7, b7 = edges[0]

    def nonsimple(w):
        w['cycles'].append([a7, b7, a7, b7])
    case('cycle_not_simple', 'appended closed walk [%d,%d,%d,%d]' % (a7, b7, a7, b7),
         nonsimple, {'cycle_not_simple'})

    # A8 a triangle (length 3) if the graph has one
    tri = None
    for (u, v) in edges:
        common = adj[u] & adj[v]
        if common:
            tri = [u, v, min(common)]
            break
    if tri:
        case('cycle_length_3', 'appended triangle %s' % tri,
             lambda w: w['cycles'].append(list(tri)), {'cycle_length'})

    # A9 a declared edge removed
    def drop_edge(w):
        w['edges'] = w['edges'][1:]
    case('edge_removed', 'first declared edge %s removed' % W['edges'][0], drop_edge,
         {'edge_mismatch'})

    # A10 a non-unit pair declared as an edge: the closest non-edge pair
    def closest_nonedge():
        best = None
        for i in range(n):
            xi, yi = fcoord(P[i], D)
            for j in range(i + 1, n):
                if (i, j) in eset:
                    continue
                xj, yj = fcoord(P[j], D)
                m = abs(math.hypot(xi - xj, yi - yj) - 1.0)
                if best is None or m < best[0]:
                    best = (m, i, j)
        return best
    m10, i10, j10 = closest_nonedge()

    def add_edge(w):
        w['edges'].append([i10, j10])
    case('edge_added', 'non-edge (%d,%d) with | |d|-1 | = %.2e declared as edge'
         % (i10, j10, m10), add_edge, {'edge_mismatch'})

    # A11 duplicate point: an off-cycle vertex copied onto a same-coloured vertex
    v11 = next(v for v in reversed(off) if v not in (v_off, fv))
    t11 = next(x for x in range(n) if x != v11 and col[x] == col[v11] and x != fv
               and x not in adj[v11])

    def dup(w):
        w['points'][v11] = list(w['points'][t11])
    case('duplicate_point', 'point %d replaced by a copy of point %d (same colour)' % (v11, t11),
         dup, {'points_not_distinct', 'edge_mismatch'})

    # A12 fixed vertex not the origin
    f12 = 0 if fv != 0 else 1

    def refix(w):
        w['fixed_vertex'] = f12
    case('fixed_not_origin', 'fixed_vertex set to %d' % f12, refix, {'fixed_not_origin'})

    # A13 non-integer coordinate
    def floaty(w):
        w['points'][0][0] = float(w['points'][0][0]) + 0.5
    case('float_coordinate', 'a coordinate replaced by a float', floaty, {'format'})

    # A14 colour out of range
    def col4(w):
        w['colouring'][0] = 4
    case('colour_out_of_range', 'a colour set to 4', col4, {'colouring_malformed'})

    log('\n== Part A: corrupted witnesses (checker run with --no-sat) ==')
    for name, desc, mutate, expected, exact, extra_check in cases:
        w = copy.deepcopy(W)
        mutate(w)
        path = os.path.join(TESTDIR, name + '.json.gz')
        with gzip.open(path, 'wt') as f:
            json.dump(w, f)
        # independent float expectation for moved points: no unit neighbours left
        if name.startswith('moved') or name.startswith('tiny'):
            v = v_on if name == 'moved_point_on_cycle' else v_off
            nb = float_unit_neighbours(w['points'], D, v)
            R.check(name + ': float view, moved vertex has no unit neighbour', nb == [],
                    '(neighbours %s)' % nb)
        jout = os.path.join(TESTDIR, name + '.out.json')
        cnf = os.path.join(TESTDIR, name + '.cnf')
        code, s, out, dt = run_checker(path, jout, cnf)
        got = set(c for c, _ in s['reasons'])
        log('\n  case %s: %s' % (name, desc))
        log('    exit code %d, %.1f s, reasons %s' % (code, dt, sorted(got)))
        for c, msg in s['reasons']:
            log('      %s: %s' % (c, msg))
        if not expected:
            R.check(name + ': accepted', code == 0 and not got and s['result'] == 'PASS')
        else:
            R.check(name + ': rejected', code == 1 and s['result'] == 'REJECT')
            if exact:
                R.check(name + ': reasons exactly %s' % sorted(expected), got == expected)
            else:
                R.check(name + ': reasons include %s' % sorted(expected), expected <= got)
        if extra_check:
            ok, what = extra_check(s)
            R.check(name + ': ' + what, ok)
        for p in (path, jout, cnf):
            if os.path.exists(p) and (p != path):
                os.remove(p)  # exact path

    # ------------------------------------------------------------------
    # Part B: encoding tests
    # ------------------------------------------------------------------
    log('\n== Part B: encoding tests ==')
    tight_given = [k for k, c in enumerate(C) if tight(col, c)]

    # B1: no cycle clauses -> SAT, model = proper colouring, fixed vertex 0.
    # Run 0 uses the formula as written.  Runs 1..9 rename the variables at
    # random, flip their signs at random and shuffle the clauses, so that the
    # solver finds different colourings; the model is mapped back.
    nv, cl, info = rc.build_cnf(n, edges, C, fv, 'arc', include_cycles=False)
    cnf1 = os.path.join(TESTDIR, 'no_cycles.cnf')
    rc.write_cnf(cnf1, nv, cl)
    log('\n  B1 formula without cycle clauses: %d vars, %d clauses, sha256 %s'
        % (nv, len(cl), rc.sha256_file(cnf1)))
    seen = set()
    min_tight = None
    for run in range(10):
        rng = random.Random(1000 + run)
        if run == 0:
            ren = {i: i for i in range(1, nv + 1)}
            sign = {i: 1 for i in range(1, nv + 1)}
            cl_run = cl
        else:
            perm = list(range(1, nv + 1))
            rng.shuffle(perm)
            ren = {i: perm[i - 1] for i in range(1, nv + 1)}
            sign = {i: rng.choice((1, -1)) for i in range(1, nv + 1)}
            cl_run = [[(1 if l > 0 else -1) * sign[abs(l)] * ren[abs(l)] for l in c] for c in cl]
            for c in cl_run:
                rng.shuffle(c)
            rng.shuffle(cl_run)
        cnfr = os.path.join(TESTDIR, 'no_cycles_run.cnf')
        rc.write_cnf(cnfr, nv, cl_run)
        code, status, model_r = kissat_solve(cnfr, run)
        os.remove(cnfr)
        # map back: old variable i is true iff new variable ren[i] has value (sign[i] == 1)
        model = {i: (model_r.get(ren[i], False) == (sign[i] == 1)) for i in range(1, nv + 1)}
        c1 = decode(model, n)
        proper = c1 is not None and all(c1[u] != c1[v] for (u, v) in edges)
        ntight = sum(1 for c in C if tight(c1, c)) if c1 else -1
        if c1:
            seen.add(tuple(c1))
            min_tight = ntight if min_tight is None else min(min_tight, ntight)
        log('    run %d: exit %s %s; decoded=%s proper=%s fixed colour=%s; listed cycles '
            'tight under this colouring: %d' % (run, code, status, c1 is not None, proper,
                                                c1[fv] if c1 else None, ntight))
        R.check('B1 run %d: SAT (exit 10)' % run, code == 10 and status == ['s SATISFIABLE'])
        R.check('B1 run %d: model is a proper 4-colouring, fixed vertex colour 0' % run,
                proper and c1[fv] == 0)
        R.check('B1 run %d: model has a tight listed cycle (as UNSAT predicts)' % run,
                ntight >= 1)
    log('    distinct colourings found: %d; fewest tight listed cycles: %s'
        % (len(seen), min_tight))
    R.check('B1: the runs found more than one colouring', len(seen) > 1)
    os.remove(cnf1)

    # B2: canonical assignment of the rotated given colouring
    rot = [(c - col[fv]) % 4 for c in col]
    for enc in ('arc', 'direct'):
        nv, cl, info = rc.build_cnf(n, edges, C, fv, enc)
        assign = {}
        for v in range(n):
            for k in range(4):
                assign[rc.xvar(v, k)] = (rot[v] == k)
        if enc == 'arc':
            # replicate the numbering: arcs in order of first appearance
            arcnum = {}
            for c in C:
                m = len(c)
                for t in range(m):
                    a = (c[t], c[(t + 1) % m])
                    if a not in arcnum:
                        arcnum[a] = 4 * n + 1 + len(arcnum)
                        assign[arcnum[a]] = (rot[a[1]] - rot[a[0]]) % 4 == 1
            R.check('B2 arc: replicated arc count equals encoder count',
                    len(arcnum) == info['arc_vars'] and nv == 4 * n + len(arcnum))
        violated = [i for i, clause in enumerate(cl)
                    if not any(assign[abs(l)] == (l > 0) for l in clause)]
        first_cycle = info['base_clauses'] + info['arc_clauses']
        non_cycle_viol = [i for i in violated if i < first_cycle]
        if enc == 'arc':
            viol_cycles = sorted(i - first_cycle for i in violated if i >= first_cycle)
        else:
            viol_cycles = sorted(set((i - first_cycle) // 4 for i in violated if i >= first_cycle))
            per = Counter((i - first_cycle) // 4 for i in violated if i >= first_cycle)
            R.check('B2 %s: each tight cycle violates exactly one of its 4 clauses' % enc,
                    set(per.values()) <= {1})
        log('\n  B2 (%s) rotated given colouring: %d clauses violated, %d of them cycle '
            'clauses; tight cycles %d' % (enc, len(violated), len(violated) - len(non_cycle_viol),
                                          len(tight_given)))
        R.check('B2 %s: all colouring/edge/arc clauses satisfied' % enc, not non_cycle_viol)
        R.check('B2 %s: violated cycle clauses == tight cycles' % enc, viol_cycles == tight_given)

    # B3: drop the tight cycles of the given colouring -> SAT, model has no tight cycle
    keep = [c for k, c in enumerate(C) if k not in set(tight_given)]
    nv, cl, info = rc.build_cnf(n, edges, keep, fv, 'arc')
    cnf3 = os.path.join(TESTDIR, 'drop_tight.cnf')
    rc.write_cnf(cnf3, nv, cl)
    code, status, model = kissat_solve(cnf3, 0)
    c3 = decode(model, n)
    proper = c3 is not None and all(c3[u] != c3[v] for (u, v) in edges)
    t3 = sum(1 for c in keep if tight(c3, c)) if c3 else -1
    log('\n  B3 %d cycles kept (%d dropped): exit %s %s; proper=%s fixed colour=%s; '
        'tight kept cycles %d' % (len(keep), len(tight_given), code, status, proper,
                                  c3[fv] if c3 else None, t3))
    R.check('B3: SAT', code == 10 and status == ['s SATISFIABLE'])
    R.check('B3: model proper, fixed vertex 0, no kept cycle tight',
            proper and c3[fv] == 0 and t3 == 0)
    os.remove(cnf3)

    log('\n%d checks, %d failures%s' % (R.count, len(R.failures),
                                        (': ' + ', '.join(R.failures)) if R.failures else ''))
    log('TESTS: %s' % ('ALL PASSED' if not R.failures else 'FAILED'))
    log.close()
    return 0 if not R.failures else 1


if __name__ == '__main__':
    sys.exit(main())
