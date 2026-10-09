#!/usr/bin/env python3
"""Independent referee check: a finite unit-distance graph H in the plane over
F = Q(sqrt3, sqrt11) with circular chromatic number 4.

Usage:
    python3 ref_check4.py WITNESS.json[.gz] OUTDIR [options]

Witness format (JSON, optionally gzipped):
    "denominator": positive integer D
    "points":  integer 8-tuples [a0,a1,a2,a3,b0,b1,b2,b3] meaning the point
               ((a0 + a1 r3 + a2 r11 + a3 r33)/D, (b0 + b1 r3 + b2 r11 + b3 r33)/D)
               where r3 = sqrt3, r11 = sqrt11, r33 = sqrt33
    "edges":   declared edges (NOT trusted; only compared at the end of task 1)
    "colouring": colours in {0,1,2,3}
    "cycles":  lists of vertex indices, each a directed cycle v0 -> v1 -> ... -> v0
    "fixed_vertex": an index

Tasks performed:
  1. recompute ALL unit-distance pairs exactly (integer arithmetic in the basis
     (1, r3, r11, r33)); check the points are distinct; compare with "edges".
  2. check the colouring is a proper 4-colouring of the recomputed graph.
  3. check each cycle is a simple closed walk along recomputed edges, length = 0 mod 4.
  4. write a CNF: "proper 4-colouring, colour(fixed_vertex) = 0, no listed cycle tight".
  5. kissat (DRAT proof) -> drat-trim (DRAT -> LRAT) -> cake_lpr (verified LRAT check).

The last line printed is "REFEREE: ACCEPTED" or "REFEREE: REJECTED <reason>"
(or "REFEREE: INCOMPLETE ..." when --no-sat is given and nothing failed).
Exit status: 0 accepted, 1 rejected, 2 incomplete.

Options:
    --no-sat               stop after task 4 (CNF written, no solver run)
    --allow-edge-mismatch  only warn (do not reject) if declared edges differ
    --keep-proofs          do not delete the DRAT/LRAT proof files at the end
    --timeout SECONDS      wall-clock limit per external tool (default 43200)
    --kissat / --drat-trim / --cake-lpr PATH   tool binaries (defaults: the environment variables
                           KISSAT, DRAT_TRIM, CAKE_LPR, else the names on the PATH)
"""

import argparse
import gzip
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import time

# tool binaries: the environment variables KISSAT, DRAT_TRIM, CAKE_LPR, else the names on the PATH
DEFAULT_KISSAT = os.environ.get("KISSAT") or shutil.which("kissat") or "kissat"
DEFAULT_DRAT_TRIM = os.environ.get("DRAT_TRIM") or shutil.which("drat-trim") or "drat-trim"
DEFAULT_CAKE_LPR = os.environ.get("CAKE_LPR") or shutil.which("cake_lpr") or "cake_lpr"

NCOL = 4


class Reject(Exception):
    """Raised as soon as the witness fails a check."""


def say(msg=""):
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Exact arithmetic in F = Q(r3, r11), elements written as integer 4-tuples in the
# basis (1, r3, r11, r33).  (The basis is linearly independent over Q because
# [F:Q] = 4: r11 is not in Q(r3), since r11 = a + b r3 with a, b rational forces
# ab = 0 and then 11 or 11/3 would be a rational square.)
# Multiplication table:  r3^2 = 3, r11^2 = 11, r33^2 = 33,
#                        r3 r11 = r33, r3 r33 = 3 r11, r11 r33 = 11 r3.
# ---------------------------------------------------------------------------
def fmul(u, w):
    u0, u1, u2, u3 = u
    w0, w1, w2, w3 = w
    return (u0 * w0 + 3 * u1 * w1 + 11 * u2 * w2 + 33 * u3 * w3,
            u0 * w1 + u1 * w0 + 11 * (u2 * w3 + u3 * w2),
            u0 * w2 + u2 * w0 + 3 * (u1 * w3 + u3 * w1),
            u0 * w3 + u3 * w0 + u1 * w2 + u2 * w1)


def fadd(u, w):
    return tuple(a + b for a, b in zip(u, w))


def fsub(u, w):
    return tuple(a - b for a, b in zip(u, w))


def dist2_numerator(p, q):
    """D^2 * |p - q|^2 as an element of F (integer 4-tuple), via the general
    multiplication table (used as an independent re-check of every edge found)."""
    dx = fsub(q[0:4], p[0:4])
    dy = fsub(q[4:8], p[4:8])
    return fadd(fmul(dx, dx), fmul(dy, dy))


def is_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


# ---------------------------------------------------------------------------
def load_witness(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return json.load(fh)


def check_structure(W):
    if not isinstance(W, dict):
        raise Reject("witness is not a JSON object")
    for key in ("denominator", "points", "colouring", "cycles", "fixed_vertex"):
        if key not in W:
            raise Reject("missing field '%s'" % key)
    # the program always reads coordinates in Q(sqrt3, sqrt11) with basis (1, r3, r11, r33);
    # if the witness declares its field / basis, make sure they agree with that reading
    import re
    for key, want in (("field", ["3", "11"]), ("basis", ["1", "3", "11", "33"])):
        if key in W:
            got = re.findall(r"\d+", str(W[key]))
            if got != want:
                raise Reject("declared %s %r is not Q(sqrt3, sqrt11) / (1, sqrt3, sqrt11, sqrt33)"
                             % (key, W[key]))
    D = W["denominator"]
    if not is_int(D) or D <= 0:
        raise Reject("denominator is not a positive integer")
    P = W["points"]
    if not isinstance(P, list) or len(P) == 0:
        raise Reject("points is not a non-empty list")
    for i, p in enumerate(P):
        if not isinstance(p, list) or len(p) != 8 or not all(is_int(x) for x in p):
            raise Reject("point %d is not an integer 8-tuple" % i)
    n = len(P)
    f = W["fixed_vertex"]
    if not is_int(f) or not (0 <= f < n):
        raise Reject("fixed_vertex is not a vertex index")
    return D, [tuple(p) for p in P], n, f


# ---------------------------------------------------------------------------
# Task 1
# ---------------------------------------------------------------------------
def unit_pairs_exact(P, D):
    """All pairs i<j with |P_i - P_j| = 1 exactly.

    With dx = (x0,x1,x2,x3), dy = (y0,..,y3) the numerators of the coordinate
    differences, D^2 |p-q|^2 = dx^2 + dy^2 has coordinates
       1   : x0^2 + 3x1^2 + 11x2^2 + 33x3^2 + (same for y)
       r3  : 2(x0x1 + 11x2x3) + 2(y0y1 + 11y2y3)
       r11 : 2(x0x2 + 3x1x3)  + 2(y0y2 + 3y1y3)
       r33 : 2(x0x3 + x1x2)   + 2(y0y3 + y1y2)
    and the distance is 1 iff this equals (D^2, 0, 0, 0)."""
    n = len(P)
    DD = D * D
    out = []
    for i in range(n):
        a0, a1, a2, a3, b0, b1, b2, b3 = P[i]
        for j in range(i + 1, n):
            q = P[j]
            x0 = q[0] - a0
            x1 = q[1] - a1
            x2 = q[2] - a2
            x3 = q[3] - a3
            y0 = q[4] - b0
            y1 = q[5] - b1
            y2 = q[6] - b2
            y3 = q[7] - b3
            if x0 * x3 + x1 * x2 + y0 * y3 + y1 * y2:
                continue
            if x0 * x1 + 11 * x2 * x3 + y0 * y1 + 11 * y2 * y3:
                continue
            if x0 * x2 + 3 * x1 * x3 + y0 * y2 + 3 * y1 * y3:
                continue
            if (x0 * x0 + 3 * x1 * x1 + 11 * x2 * x2 + 33 * x3 * x3 +
                    y0 * y0 + 3 * y1 * y1 + 11 * y2 * y2 + 33 * y3 * y3) != DD:
                continue
            out.append((i, j))
    return out


def float_crosscheck(P, D, exact_set):
    """Independent floating-point view of all pairs (numpy, double precision):
    every pair with |d^2 - 1| < 1e-6 must be an exact unit pair and vice versa.
    Returns (n_candidates, n_mismatch, min |d^2-1| over non-unit pairs) or None."""
    try:
        import numpy as np
    except ImportError:
        return None
    r = np.array([1.0, 3.0 ** 0.5, 11.0 ** 0.5, 33.0 ** 0.5])
    A = np.array(P, dtype=np.float64)
    X = A[:, 0:4] @ r / D
    Y = A[:, 4:8] @ r / D
    n = len(P)
    cand = set()
    gap = float("inf")
    for i in range(n - 1):
        d2 = (X[i + 1:] - X[i]) ** 2 + (Y[i + 1:] - Y[i]) ** 2
        err = np.abs(d2 - 1.0)
        for k in np.nonzero(err < 1e-6)[0]:
            cand.add((i, i + 1 + int(k)))
        mask = np.ones(len(err), dtype=bool)
        for k in np.nonzero(err < 1e-6)[0]:
            mask[k] = False
        if mask.any():
            gap = min(gap, float(err[mask].min()))
    mism = len(cand ^ exact_set)
    return len(cand), mism, gap


def task1(W, D, P, n, args, summary):
    say("== Task 1: points and unit-distance graph ==")
    say("points n = %d, denominator D = %d, max |integer coordinate| = %d"
        % (n, D, max(abs(x) for p in P for x in p)))
    # distinctness: with a common denominator and a Q-basis, two points are equal
    # iff their integer 8-tuples are equal
    if len(set(P)) != n:
        seen = {}
        for i, p in enumerate(P):
            if p in seen:
                raise Reject("points %d and %d coincide" % (seen[p], i))
            seen[p] = i
    say("all %d points are distinct" % n)
    t0 = time.time()
    E = unit_pairs_exact(P, D)
    say("exact all-pairs scan: %d pairs examined, %d unit-distance pairs (%.1f s)"
        % (n * (n - 1) // 2, len(E), time.time() - t0))
    # re-check every edge with the general multiplication table
    for (i, j) in E:
        if dist2_numerator(P[i], P[j]) != (D * D, 0, 0, 0):
            raise Reject("internal inconsistency on pair %d %d" % (i, j))
    Eset = set(E)
    fc = float_crosscheck(P, D, Eset)
    if fc is None:
        say("float cross-check skipped (numpy not available)")
    else:
        say("float cross-check: %d pairs with |d^2-1| < 1e-6, %d disagreements with "
            "the exact set; min |d^2-1| over the other pairs = %.3g" % fc)
        if fc[1]:
            raise Reject("exact and floating-point unit-pair sets disagree")
    deg = [0] * n
    for (i, j) in E:
        deg[i] += 1
        deg[j] += 1
    say("degrees: min %d, max %d, isolated vertices %d"
        % (min(deg), max(deg), sum(1 for d in deg if d == 0)))
    summary["n_points"] = n
    summary["n_edges_recomputed"] = len(E)
    summary["float_crosscheck"] = fc
    return E, Eset


def compare_declared_edges(W, n, Eset, args, summary):
    say("== Task 1 (end): comparison with the declared 'edges' field ==")
    if "edges" not in W:
        say("no declared 'edges' field; nothing to compare")
        summary["edges_compare"] = "absent"
        return
    decl = set()
    bad = 0
    dup = 0
    for e in W["edges"]:
        if (not isinstance(e, list) or len(e) != 2 or not all(is_int(x) for x in e)
                or not (0 <= e[0] < n and 0 <= e[1] < n) or e[0] == e[1]):
            bad += 1
            continue
        key = (min(e), max(e))
        if key in decl:
            dup += 1
        decl.add(key)
    missing = sorted(Eset - decl)   # unit pairs not declared
    extra = sorted(decl - Eset)     # declared pairs that are not unit pairs
    say("declared entries %d (malformed %d, duplicates %d, distinct %d); recomputed %d"
        % (len(W["edges"]), bad, dup, len(decl), len(Eset)))
    say("unit pairs missing from declared list: %d %s" % (len(missing), missing[:10]))
    say("declared pairs that are NOT at unit distance: %d %s" % (len(extra), extra[:10]))
    summary["edges_compare"] = {"declared": len(W["edges"]), "malformed": bad,
                                "duplicates": dup, "missing": len(missing),
                                "extra": len(extra)}
    if bad or dup or missing or extra:
        if args.allow_edge_mismatch:
            say("WARNING: declared edges differ from the recomputed graph "
                "(ignored because of --allow-edge-mismatch; all checks use the recomputed graph)")
        else:
            raise Reject("declared edges differ from recomputed unit-distance graph "
                         "(malformed %d, duplicates %d, missing %d, non-unit %d)"
                         % (bad, dup, len(missing), len(extra)))
    else:
        say("declared edge set == recomputed edge set")


# ---------------------------------------------------------------------------
# Task 2
# ---------------------------------------------------------------------------
def task2(W, n, E, summary):
    say("== Task 2: the colouring ==")
    c = W["colouring"]
    if not isinstance(c, list) or len(c) != n:
        raise Reject("colouring does not have one entry per point")
    for v, k in enumerate(c):
        if not is_int(k) or not (0 <= k < NCOL):
            raise Reject("colour of vertex %d is not in {0,1,2,3}" % v)
    for (u, v) in E:
        if c[u] == c[v]:
            raise Reject("colouring not proper: edge %d-%d both coloured %d" % (u, v, c[u]))
    counts = [c.count(k) for k in range(NCOL)]
    say("proper 4-colouring of the recomputed graph; colour class sizes %s" % counts)
    summary["colour_classes"] = counts
    return c


# ---------------------------------------------------------------------------
# Task 3
# ---------------------------------------------------------------------------
def task3(W, n, Eset, summary):
    say("== Task 3: the listed cycles ==")
    cyc = W["cycles"]
    if not isinstance(cyc, list) or len(cyc) == 0:
        raise Reject("no cycles listed")
    lengths = {}
    for idx, C in enumerate(cyc):
        if not isinstance(C, list) or not all(is_int(x) for x in C):
            raise Reject("cycle %d is not a list of integers" % idx)
        m = len(C)
        if m < 3 or m % 4 != 0:
            raise Reject("cycle %d has length %d (not a positive multiple of 4)" % (idx, m))
        if not all(0 <= x < n for x in C):
            raise Reject("cycle %d contains an index out of range" % idx)
        if len(set(C)) != m:
            raise Reject("cycle %d repeats a vertex (not simple)" % idx)
        for i in range(m):
            a, b = C[i], C[(i + 1) % m]
            if (min(a, b), max(a, b)) not in Eset:
                raise Reject("cycle %d uses non-edge %d -> %d" % (idx, a, b))
        lengths[m] = lengths.get(m, 0) + 1
    arcs = set()
    for C in cyc:
        for i in range(len(C)):
            arcs.add((C[i], C[(i + 1) % len(C)]))
    verts = set(x for C in cyc for x in C)
    say("%d cycles, all simple closed walks along recomputed edges; lengths %s"
        % (len(cyc), dict(sorted(lengths.items()))))
    say("distinct directed arcs used: %d; vertices touched: %d" % (len(arcs), len(verts)))
    dup = len(cyc) - len(set(tuple(C) for C in cyc))
    if dup:
        say("note: %d listed cycles are exact duplicates of others (harmless)" % dup)
    summary["n_cycles"] = len(cyc)
    summary["cycle_lengths"] = {str(k): v for k, v in sorted(lengths.items())}
    summary["n_arcs"] = len(arcs)
    return cyc


def tight_cycles(c, cycles):
    """Indices of listed cycles all of whose arcs a->b satisfy c(b)-c(a) = 1 mod 4."""
    res = []
    for idx, C in enumerate(cycles):
        m = len(C)
        if all((c[C[(i + 1) % m]] - c[C[i]]) % NCOL == 1 for i in range(m)):
            res.append(idx)
    return res


# ---------------------------------------------------------------------------
# Task 4: CNF
#   x(v,k) = k*n + v + 1          (v in [0,n), k in {0,1,2,3}): vertex v has colour k
#   t(a,b) = 4n + 1 + index       one variable per directed arc a->b used by a cycle
# Clauses
#   (A) x(v,0) | x(v,1) | x(v,2) | x(v,3)                      each v
#   (B) -x(v,k) | -x(v,l)                                      each v, k<l
#   (C) -x(u,k) | -x(v,k)                                      each edge uv, each k
#   (D) x(f,0)                                                 fixed vertex
#   (E) -x(a,k) | -x(b,k+1 mod 4) | t(a,b)                     each arc, each k
#   (F) OR_{arcs a->b of the cycle} -t(a,b)                    each listed cycle
# ---------------------------------------------------------------------------
def build_cnf(n, E, cycles, f, drop_cycle_clauses=False):
    def x(v, k):
        return k * n + v + 1

    arc_index = {}
    for C in cycles:
        m = len(C)
        for i in range(m):
            a = (C[i], C[(i + 1) % m])
            if a not in arc_index:
                arc_index[a] = len(arc_index)

    def t(a, b):
        return NCOL * n + 1 + arc_index[(a, b)]

    nvars = NCOL * n + len(arc_index)
    cls = []
    for v in range(n):
        cls.append([x(v, k) for k in range(NCOL)])
    for v in range(n):
        for k in range(NCOL):
            for l in range(k + 1, NCOL):
                cls.append([-x(v, k), -x(v, l)])
    for (u, v) in E:
        for k in range(NCOL):
            cls.append([-x(u, k), -x(v, k)])
    cls.append([x(f, 0)])
    for (a, b) in arc_index:
        for k in range(NCOL):
            cls.append([-x(a, k), -x(b, (k + 1) % NCOL), t(a, b)])
    if not drop_cycle_clauses:
        for C in cycles:
            m = len(C)
            cls.append([-t(C[i], C[(i + 1) % m]) for i in range(m)])
    meta = {"x": x, "t": t, "arc_index": arc_index, "nvars": nvars}
    return nvars, cls, meta


def write_cnf(path, nvars, cls, comment):
    with open(path, "w") as fh:
        for line in comment:
            fh.write("c %s\n" % line)
        fh.write("p cnf %d %d\n" % (nvars, len(cls)))
        for c in cls:
            fh.write(" ".join(map(str, c)) + " 0\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Task 5: external tools
# ---------------------------------------------------------------------------
def run_tool(cmd, logpath, timeout):
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.time()
    with open(logpath, "w") as fh:
        try:
            p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
            rc = p.returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    wall = time.time() - t0
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (r1.ru_utime - r0.ru_utime) + (r1.ru_stime - r0.ru_stime)
    with open(logpath, errors="replace") as fh:
        lines = fh.read().splitlines()
    slines = [l.strip() for l in lines if l.startswith("s ")]
    return rc, slines, wall, cpu, lines


def remove_exact(path):
    if os.path.isfile(path):
        os.remove(path)


def task5(cnf, outdir, args, summary):
    say("== Task 5: kissat -> drat-trim (LRAT) -> cake_lpr ==")
    drat = os.path.join(outdir, "ref4_proof.drat")
    lrat = os.path.join(outdir, "ref4_proof.lrat")
    remove_exact(drat)
    remove_exact(lrat)
    h0 = sha256(cnf)
    try:
        rc, s, wall, cpu, _ = run_tool([args.kissat, cnf, drat],
                                       os.path.join(outdir, "kissat.log"), args.timeout)
        say("kissat: exit %s, %s, wall %.1f s, cpu %.1f s, proof %d bytes"
            % (rc, s, wall, cpu, os.path.getsize(drat) if os.path.exists(drat) else -1))
        summary["kissat"] = {"exit": rc, "s": s, "wall": wall, "cpu": cpu}
        if s != ["s UNSATISFIABLE"]:
            raise Reject("kissat did not report UNSATISFIABLE (%s)" % s)
        rc, s, wall, cpu, lines = run_tool([args.drat_trim, cnf, drat, "-L", lrat],
                                           os.path.join(outdir, "drat_trim.log"), args.timeout)
        say("drat-trim: exit %s, %s, wall %.1f s, cpu %.1f s, LRAT %d bytes"
            % (rc, s, wall, cpu, os.path.getsize(lrat) if os.path.exists(lrat) else -1))
        summary["drat_trim"] = {"exit": rc, "s": s, "wall": wall, "cpu": cpu}
        if "s VERIFIED" not in s:
            raise Reject("drat-trim did not verify the proof (%s)" % s)
        rc, s, wall, cpu, lines = run_tool([args.cake_lpr, cnf, lrat],
                                           os.path.join(outdir, "cake_lpr.log"), args.timeout)
        say("cake_lpr: exit %s, %s, wall %.1f s, cpu %.1f s" % (rc, s, wall, cpu))
        summary["cake_lpr"] = {"exit": rc, "s": s, "wall": wall, "cpu": cpu}
        if s != ["s VERIFIED UNSAT"] or rc != 0:
            raise Reject("cake_lpr did not print 's VERIFIED UNSAT' (%s, exit %s)" % (s, rc))
    finally:
        if not args.keep_proofs:
            remove_exact(drat)
            remove_exact(lrat)
    h1 = sha256(cnf)
    if h0 != h1:
        raise Reject("CNF file changed while the tools ran")
    say("CNF sha256 %s unchanged during the run" % h0)


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("witness")
    ap.add_argument("outdir")
    ap.add_argument("--no-sat", action="store_true")
    ap.add_argument("--allow-edge-mismatch", action="store_true")
    ap.add_argument("--keep-proofs", action="store_true")
    ap.add_argument("--timeout", type=float, default=43200.0)
    ap.add_argument("--kissat", default=DEFAULT_KISSAT)
    ap.add_argument("--drat-trim", dest="drat_trim", default=DEFAULT_DRAT_TRIM)
    ap.add_argument("--cake-lpr", dest="cake_lpr", default=DEFAULT_CAKE_LPR)
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    outdir = os.path.abspath(args.outdir)
    summary = {"witness": os.path.abspath(args.witness)}
    try:
        try:
            W = load_witness(args.witness)
        except (OSError, ValueError) as e:
            raise Reject("cannot read witness: %s" % e)
        D, P, n, f = check_structure(W)
        E, Eset = task1(W, D, P, n, args, summary)
        c = task2(W, n, E, summary)
        cycles = task3(W, n, Eset, summary)
        compare_declared_edges(W, n, Eset, args, summary)
        tc = tight_cycles(c, cycles)
        say("(consistency) listed cycles tight under the given colouring: %d" % len(tc))
        say("== Task 4: CNF ==")
        nvars, cls, meta = build_cnf(n, E, cycles, f)
        cnf = os.path.join(outdir, "ref4_main.cnf")
        write_cnf(cnf, nvars, cls, [
            "referee CNF: proper 4-colouring, colour(fixed vertex %d) = 0, no listed cycle tight" % f,
            "x(v,k) = k*%d + v + 1 ; t(arc) = %d + 1 + arc index" % (n, NCOL * n),
            "witness %s" % os.path.basename(args.witness)])
        say("CNF %s: %d variables (%d colour + %d tight-arc), %d clauses, sha256 %s"
            % (cnf, nvars, NCOL * n, len(meta["arc_index"]), len(cls), sha256(cnf)))
        summary["cnf"] = {"path": cnf, "vars": nvars, "clauses": len(cls), "sha256": sha256(cnf)}
        if args.no_sat:
            verdict = "REFEREE: INCOMPLETE (tasks 1-4 passed, SAT stage skipped by --no-sat)"
            code = 2
        else:
            task5(cnf, outdir, args, summary)
            verdict = "REFEREE: ACCEPTED"
            code = 0
    except Reject as e:
        verdict = "REFEREE: REJECTED %s" % e
        code = 1
    summary["verdict"] = verdict
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=str)
    say(verdict)
    return code


if __name__ == "__main__":
    sys.exit(main())
