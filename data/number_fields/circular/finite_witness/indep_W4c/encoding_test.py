#!/usr/bin/env python3
"""Semantic test of the CNF FILE actually given to the solver/checkers.

For many colourings c: V -> {0,1,2,3} we compare, on the assignment
x(v,k) = [c(v) == k]  (x(v,k) = k*n + v + 1):

   CNF_viol(c)    = number of clauses of the CNF file falsified, and
   DIRECT_viol(c) = (#monochromatic listed edges) + [c(fixed) != 0]
                    + (#listed cycles tight in their listed direction,
                       evaluated literally: c(v_{t+1}) - c(v_t) = 1 mod 4, all t)
                    (the last term is omitted with --no-cycles).

They must be EQUAL for every colouring tried (each monochromatic edge
falsifies exactly one edge clause, each tight cycle exactly one of its four
shift clauses).  In particular CNF(c) is true iff DIRECT(c) is true.
Since the ALO/AMO clauses force every model to be of the form x_c, this
pins down the semantics of the file.  Sample:
  * the stored witness colouring under all 24 colour permutations, each also
    renormalised so that c(fixed) = 0 by a cyclic shift (tightness-preserving);
  * a colouring decoded from a model file (optional), same treatment;
  * random single-vertex recolourings of those (mostly improper);
  * purely random colourings.
Also checks that an assignment with two colours on a vertex / no colour on a
vertex is rejected by the CNF.

Usage: python3 encoding_test.py WITNESS.json.gz CNF [--model KISSAT_OUT] [--no-cycles]
"""
import gzip
import itertools
import json
import random
import sys

random.seed(77)


def read_cnf(path):
    clauses, cur = [], []
    nv = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("p"):
                nv = int(line.split()[2])
                continue
            if line.startswith("c"):
                continue
            for tok in line.split():
                l = int(tok)
                if l == 0:
                    clauses.append(cur)
                    cur = []
                else:
                    cur.append(l)
    assert not cur
    return nv, clauses


def main():
    args = sys.argv[1:]
    with gzip.open(args[0], "rt") as fh:
        d = json.load(fh)
    nv, clauses = read_cnf(args[1])
    use_cycles = "--no-cycles" not in args
    n = len(d["points"])
    assert nv == 4 * n
    E = [tuple(e) for e in d["edges"]]
    cycles = d["cycles"] if use_cycles else []
    fv = d["fixed_vertex"]

    def n_tight(c):
        k = 0
        for cyc in cycles:
            m = len(cyc)
            if all((c[cyc[(t + 1) % m]] - c[cyc[t]]) % 4 == 1 for t in range(m)):
                k += 1
        return k

    def direct_viol(c):
        mono = sum(1 for (u, v) in E if c[u] == c[v])
        return mono + (1 if c[fv] != 0 else 0) + n_tight(c)

    def cnf_viol_assign(val):
        return sum(1 for cl in clauses
                   if not any((val[abs(l)] if l > 0 else not val[abs(l)]) for l in cl))

    def cnf_viol(c):
        val = {}
        for v in range(n):
            for k in range(4):
                val[k * n + v + 1] = (c[v] == k)
        return cnf_viol_assign(val)

    base_cols = [("stored", d["colouring"])]
    if "--model" in args:
        lits = []
        with open(args[args.index("--model") + 1]) as fh:
            for line in fh:
                if line.startswith("v "):
                    lits.extend(int(t) for t in line[2:].split())
        pos = set(l for l in lits if l > 0)
        mc = []
        for v in range(n):
            ks = [k for k in range(4) if k * n + v + 1 in pos]
            assert len(ks) == 1
            mc.append(ks[0])
        base_cols.append(("model", mc))

    stats = {"tested": 0, "agree": 0, "cnf_true": 0, "direct_true": 0,
             "proper_fixed0": 0, "min_tight_among_proper_fixed0": None}
    disagreements = []
    viols = []

    def test(c, tag):
        va, vb = cnf_viol(c), direct_viol(c)
        stats["tested"] += 1
        viols.append(vb)
        if va == vb:
            stats["agree"] += 1
        else:
            disagreements.append((tag, va, vb))
        stats["cnf_true"] += (va == 0)
        stats["direct_true"] += (vb == 0)
        if c[fv] == 0 and all(c[u] != c[v] for (u, v) in E):
            stats["proper_fixed0"] += 1
            t = n_tight(c)
            m = stats["min_tight_among_proper_fixed0"]
            stats["min_tight_among_proper_fixed0"] = t if m is None else min(m, t)

    for name, col in base_cols:
        for perm in itertools.permutations(range(4)):
            c = [perm[x] for x in col]
            test(c, (name, perm, "raw"))
            sh = [(x - c[fv]) % 4 for x in c]
            test(sh, (name, perm, "shifted"))
            for _ in range(5):
                c2 = list(sh)
                v = random.randrange(n)
                c2[v] = random.randrange(4)
                test(c2, (name, perm, "recolour"))
    for _ in range(50):
        c = [random.randrange(4) for _ in range(n)]
        c[fv] = 0
        test(c, ("random",))
    # ALO / AMO behaviour on non-colouring assignments
    col = [(x - d["colouring"][fv]) % 4 for x in d["colouring"]]
    val = {k * n + v + 1: (col[v] == k) for v in range(n) for k in range(4)}
    v0 = (fv + 1) % n
    val2 = dict(val)
    val2[((col[v0] + 1) % 4) * n + v0 + 1] = True     # two colours on v0
    val3 = dict(val)
    val3[col[v0] * n + v0 + 1] = False                 # no colour on v0
    def falsified_set(vl):
        return set(i for i, cl in enumerate(clauses)
                   if not any((vl[abs(l)] if l > 0 else not vl[abs(l)]) for l in cl))
    base_set = falsified_set(val)
    # the corrupted assignment must falsify some clause that the colouring itself satisfies
    amo_rejected = bool(falsified_set(val2) - base_set)
    alo_rejected = bool(falsified_set(val3) - base_set)
    print("CNF file:", args[1], "| cycles in DIRECT:", use_cycles)
    print("colourings tested: %d; CNF_viol == DIRECT_viol on %d; disagreements %d"
          % (stats["tested"], stats["agree"], len(disagreements)))
    print("CNF-true: %d, DIRECT-true: %d; violation counts range %d..%d"
          % (stats["cnf_true"], stats["direct_true"], min(viols), max(viols)))
    print("proper colourings with c(fixed)=0 among them: %d; min #tight listed cycles among those: %s"
          % (stats["proper_fixed0"], stats["min_tight_among_proper_fixed0"]))
    print("two-colours-on-a-vertex assignment falsifies a new clause:", amo_rejected)
    print("no-colour-on-a-vertex assignment falsifies a new clause:", alo_rejected)
    if disagreements:
        print("first disagreements:", disagreements[:5])
    ok = not disagreements and amo_rejected and alo_rejected
    print("ENCODING TEST:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
