#!/usr/bin/env python3
"""Sanity tests for ref_check4.py.

Usage:  python3 ref_tests4.py WITNESS.json[.gz] OUTDIR

 (a1) one coordinate of one point (a vertex on a listed cycle) changed      -> must REJECT
 (a2) same, with the declared-edge comparison switched off                  -> must REJECT
 (a3) one coordinate of a point lying on no listed cycle changed            -> must REJECT
      (only the comparison with the declared edge list can see this one)
 (b)  colouring with two adjacent vertices given the same colour            -> must REJECT
 (c)  a listed cycle rerouted through a non-edge                            -> must REJECT
 (d)  CNF without the cycle clauses                                         -> must be SAT,
      and the decoded model must be a proper 4-colouring with colour(fixed) = 0
 (e)  encoding self-check: the assignment induced by the given colouring (rotated so that
      colour(fixed) = 0) violates exactly the cycle clauses of the cycles that are tight
      under that colouring, and nothing else
Prints "TESTS: ALL PASSED" or "TESTS: FAILED ..." as the last line.
"""

import gzip
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_check4 as R  # noqa: E402

CHECKER = os.path.join(HERE, "ref_check4.py")


def save(W, path):
    with gzip.open(path, "wt") as fh:
        json.dump(W, fh)


def run_checker(wpath, outdir, extra=()):
    cmd = [sys.executable, CHECKER, wpath, outdir, "--no-sat"] + list(extra)
    p = subprocess.run(cmd, capture_output=True, text=True)
    lines = p.stdout.strip().splitlines()
    return p.returncode, (lines[-1] if lines else "(no output)")


def main():
    wpath, outdir = sys.argv[1], os.path.abspath(sys.argv[2])
    os.makedirs(outdir, exist_ok=True)
    W = R.load_witness(wpath)
    D, P, n, f = R.check_structure(W)
    E = R.unit_pairs_exact(P, D)
    Eset = set(E)
    adj = [set() for _ in range(n)]
    for (u, v) in E:
        adj[u].add(v)
        adj[v].add(u)
    cycles = W["cycles"]
    failures = []

    def expect_reject(name, Wmod, extra=(), must_contain=None):
        p = os.path.join(outdir, "test_%s.json.gz" % name)
        save(Wmod, p)
        rc, last = run_checker(p, os.path.join(outdir, "out_%s" % name), extra)
        ok = last.startswith("REFEREE: REJECTED") and rc == 1
        if ok and must_contain is not None and must_contain not in last:
            ok = False
        print("(%s) %s  [exit %d]  -> %s" % (name, last, rc, "PASS" if ok else "FAIL"), flush=True)
        if not ok:
            failures.append(name)
        os.remove(p)

    # (a) change one coordinate of one point lying on a listed cycle
    v = cycles[0][1] if cycles[0][1] != f else cycles[0][2]
    Wa = json.loads(json.dumps(W))
    Wa["points"][v][0] += 1
    print("(a) vertex %d: coordinate a0 %d -> %d" % (v, W["points"][v][0], Wa["points"][v][0]))
    expect_reject("a1", Wa)
    expect_reject("a2", Wa, ["--allow-edge-mismatch"], must_contain="non-edge")
    # (a3) same kind of change on a vertex that lies on NO listed cycle: only the
    # comparison with the declared edge list can notice it (strict mode)
    on_cycle = set(x for C in cycles for x in C)
    v3 = next(x for x in range(n) if x not in on_cycle and x != f and adj[x])
    Wa3 = json.loads(json.dumps(W))
    Wa3["points"][v3][0] += 1
    print("(a3) vertex %d (on no listed cycle, degree %d): coordinate a0 %d -> %d"
          % (v3, len(adj[v3]), W["points"][v3][0], Wa3["points"][v3][0]))
    expect_reject("a3", Wa3, must_contain="declared edges differ")

    # (b) two adjacent vertices with the same colour
    u, w = E[0]
    Wb = json.loads(json.dumps(W))
    Wb["colouring"][w] = Wb["colouring"][u]
    print("(b) edge %d-%d: colour of %d set to %d" % (u, w, w, Wb["colouring"][w]))
    expect_reject("b", Wb, must_contain="not proper")

    # (c) a listed cycle rerouted through a non-edge
    C = list(cycles[0])
    z = next(z for z in range(n) if z not in adj[C[0]] and z not in C and z != C[0])
    Wc = json.loads(json.dumps(W))
    Wc["cycles"][0] = [C[0], z] + C[2:]
    print("(c) cycle 0 %s -> %s (%d-%d is not an edge)" % (C, Wc["cycles"][0], C[0], z))
    expect_reject("c", Wc, must_contain="non-edge")

    # (d) drop all cycle clauses -> SAT, decode and check the model
    nvars, cls, meta = R.build_cnf(n, E, cycles, f, drop_cycle_clauses=True)
    cnf = os.path.join(outdir, "nocycles.cnf")
    R.write_cnf(cnf, nvars, cls, ["referee CNF without cycle clauses"])
    rc, s, wall, cpu, lines = R.run_tool([R.DEFAULT_KISSAT, cnf],
                                         os.path.join(outdir, "kissat_nocycles.log"), 3600)
    print("(d) CNF without cycle clauses: %d vars, %d clauses; kissat exit %s %s (wall %.1f s, cpu %.1f s)"
          % (nvars, len(cls), rc, s, wall, cpu))
    okd = s == ["s SATISFIABLE"] and rc == 10
    if okd:
        val = {}
        for l in lines:
            if l.startswith("v "):
                for tok in l[2:].split():
                    lit = int(tok)
                    if lit:
                        val[abs(lit)] = lit > 0
        # every clause satisfied by the model?
        okd = all(any(val.get(abs(l), False) == (l > 0) for l in c) for c in cls)
        col = []
        for x in range(n):
            ks = [k for k in range(R.NCOL) if val.get(meta["x"](x, k), False)]
            if len(ks) != 1:
                okd = False
                ks = [None]
            col.append(ks[0])
        proper = okd and all(col[a] != col[b] for (a, b) in E)
        tc = R.tight_cycles(col, cycles) if proper else []
        print("(d) model satisfies all clauses: %s; decoded colouring proper: %s; colour(fixed) = %s; "
              "listed cycles tight under it: %d" % (okd, proper, col[f], len(tc)))
        okd = okd and proper and col[f] == 0
    print("(d) -> %s" % ("PASS" if okd else "FAIL"), flush=True)
    if not okd:
        failures.append("d")
    os.remove(cnf)

    # (e) encoding self-check with the given colouring
    c0 = W["colouring"]
    cr = [(k - c0[f]) % R.NCOL for k in c0]
    nvars, cls, meta = R.build_cnf(n, E, cycles, f)
    val = {}
    for x in range(n):
        for k in range(R.NCOL):
            val[meta["x"](x, k)] = (cr[x] == k)
    for (a, b), _ in meta["arc_index"].items():
        val[meta["t"](a, b)] = ((cr[b] - cr[a]) % R.NCOL == 1)
    violated = [i for i, c in enumerate(cls) if not any(val[abs(l)] == (l > 0) for l in c)]
    first_cycle_clause = len(cls) - len(cycles)
    viol_cycles = sorted(i - first_cycle_clause for i in violated)
    tc = R.tight_cycles(cr, cycles)
    tc_orig = R.tight_cycles(c0, cycles)
    # rotation invariance, checked for every shift s: c+s is proper iff c is, same tight cycles
    for s in range(1, R.NCOL):
        cs = [(k + s) % R.NCOL for k in c0]
        if R.tight_cycles(cs, cycles) != tc_orig or any(cs[a] == cs[b] for (a, b) in E):
            tc_orig = None
    oke = all(i >= first_cycle_clause for i in violated) and viol_cycles == tc and tc == tc_orig
    print("(e) induced assignment of the given colouring (rotated by %d): %d violated clauses, all of "
          "them cycle clauses: %s; they are exactly the %d tight listed cycles: %s; every rotation c+s "
          "keeps properness and the tight set: %s -> %s"
          % (-c0[f] % R.NCOL, len(violated), all(i >= first_cycle_clause for i in violated),
             len(tc), viol_cycles == tc, tc == tc_orig, "PASS" if oke else "FAIL"), flush=True)
    if not oke:
        failures.append("e")

    print("TESTS: ALL PASSED" if not failures else "TESTS: FAILED %s" % failures)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
