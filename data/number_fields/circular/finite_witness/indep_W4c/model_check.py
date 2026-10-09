#!/usr/bin/env python3
"""Decode / evaluate colourings against the witness.

Usage:
  python3 model_check.py WITNESS.json.gz --model KISSAT_OUTPUT_FILE
      decode the SAT model of the referee CNF (x(v,c) = c*n+v+1), check it
      is a proper 4-colouring with c(fixed)=0, and count its tight listed
      cycles (it must have >= 1 if the UNSAT claim holds).
  python3 model_check.py WITNESS.json.gz --witness-colouring
      count tight listed cycles of the colouring stored in the witness
      (again it must have >= 1), and evaluate every clause of the referee
      encoding on the shifted colouring c - c(fixed).
"""
import gzip
import json
import sys


def tight_cycles(cycles, col):
    out = []
    for k, cyc in enumerate(cycles):
        m = len(cyc)
        if all((col[cyc[(t + 1) % m]] - col[cyc[t]]) % 4 == 1 for t in range(m)):
            out.append(k)
    return out


def main():
    args = sys.argv[1:]
    with gzip.open(args[0], "rt") as fh:
        d = json.load(fh)
    n = len(d["points"])
    E = [tuple(e) for e in d["edges"]]
    fv = d["fixed_vertex"]
    if "--model" in args:
        path = args[args.index("--model") + 1]
        status = None
        lits = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("s "):
                    status = line.strip()
                elif line.startswith("v "):
                    lits.extend(int(t) for t in line[2:].split())
        print("solver status line:", status)
        assert status == "s SATISFIABLE"
        val = {}
        for l in lits:
            if l != 0:
                val[abs(l)] = l > 0
        col = []
        for v in range(n):
            cs = [c for c in range(4) if val.get(c * n + v + 1, False)]
            assert len(cs) == 1, ("vertex", v, cs)
            col.append(cs[0])
        bad = sum(1 for (u, v) in E if col[u] == col[v])
        print("decoded colouring: monochromatic edges =", bad, " c(fixed) =", col[fv])
        tc = tight_cycles(d["cycles"], col)
        print("tight listed cycles in this model:", len(tc), "(first few:", tc[:5], ")")
        ok = bad == 0 and col[fv] == 0 and len(tc) >= 1
        print("MODEL CHECK:", "PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)
    if "--witness-colouring" in args:
        col = d["colouring"]
        tc = tight_cycles(d["cycles"], col)
        sh = [(c - col[fv]) % 4 for c in col]
        tc2 = tight_cycles(d["cycles"], sh)
        print("stored colouring: c(fixed) =", col[fv], "; tight listed cycles:", len(tc),
              "; after shift to c(fixed)=0:", len(tc2), "; first few:", tc[:5])
        # evaluate the referee encoding clauses directly (semantic check of the encoding)
        violated = {"alo": 0, "amo": 0, "edge": 0, "unit": 0, "cycle": 0}
        X = lambda v, c: sh[v] == c
        for v in range(n):
            if not any(X(v, c) for c in range(4)):
                violated["alo"] += 1
            for c in range(4):
                for c2 in range(c + 1, 4):
                    if X(v, c) and X(v, c2):
                        violated["amo"] += 1
        for (u, v) in E:
            for c in range(4):
                if X(u, c) and X(v, c):
                    violated["edge"] += 1
        if not X(fv, 0):
            violated["unit"] += 1
        viol_cyc = set()
        for k, cyc in enumerate(d["cycles"]):
            m = len(cyc)
            for s in range(4):
                if all(X(cyc[t], (s + t) % 4) for t in range(m)):
                    violated["cycle"] += 1
                    viol_cyc.add(k)
        print("referee clauses violated by shifted stored colouring:", violated)
        print("cycles with a violated cycle clause == tight cycles:", viol_cyc == set(tc2))
        ok = len(tc) >= 1 and violated["alo"] == violated["amo"] == violated["edge"] == violated["unit"] == 0 \
            and viol_cyc == set(tc2) and tc == tc2
        print("WITNESS COLOURING CHECK:", "PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
