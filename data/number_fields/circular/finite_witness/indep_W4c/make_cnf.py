#!/usr/bin/env python3
"""Referee's own CNF encoding of

   "c is a proper 4-colouring of H, c(fixed_vertex) = 0, and no listed
    cycle is tight in its listed direction".

Usage: python3 make_cnf.py WITNESS.json.gz OUT.cnf [--no-cycles] [--map OUT.map.json]

Variables (colour-major numbering, chosen on purpose to differ from a
vertex-major layout):  x(v, c) = c*n + v + 1  means "vertex v has colour c",
v in 0..n-1, c in 0..3.

Clauses
  (1) for every vertex v:           x(v,0) | x(v,1) | x(v,2) | x(v,3)
  (2) for every vertex v, c < c':   -x(v,c) | -x(v,c')
  (3) for every edge uv, colour c:  -x(u,c) | -x(v,c)
  (4) unit clause                   x(fixed_vertex, 0)
  (5) for every listed cycle (v_0..v_{m-1}) and every shift s in 0..3:
          OR_{t=0}^{m-1}  -x(v_t, (s + t) mod 4)
      A colouring is tight on the cycle in the listed direction iff
      c(v_t) = c(v_0) + t (mod 4) for all t and m = 0 (mod 4); with
      s = c(v_0) the clause for shift s is exactly the negation of that.
      A cycle with m != 0 mod 4 can never be tight, so it gets no clause
      (none occur here; the geometry checker verifies m = 0 mod 4).

The edge set used is the listed edge set; this script re-derives the
Cayley pairs from points and generators and refuses to continue if the two
sets differ.
"""
import gzip
import json
import sys


def main():
    args = sys.argv[1:]
    src, out = args[0], args[1]
    no_cycles = "--no-cycles" in args
    with gzip.open(src, "rt") as fh:
        d = json.load(fh)
    P = [tuple(p) for p in d["points"]]
    n = len(P)
    assert len(set(P)) == n
    gens = set()
    for g in d["generators"]:
        gens.add(tuple(g))
        gens.add(tuple(-x for x in g))
    idx = {p: i for i, p in enumerate(P)}
    cay = set()
    for i, p in enumerate(P):
        for g in gens:
            j = idx.get(tuple(a + b for a, b in zip(p, g)))
            if j is not None:
                cay.add((min(i, j), max(i, j)))
    E = set((min(e), max(e)) for e in d["edges"])
    assert E == cay, "listed edges differ from Cayley pairs"
    assert len(E) == len(d["edges"])

    def x(v, c):
        return c * n + v + 1

    clauses = []
    for v in range(n):
        clauses.append([x(v, 0), x(v, 1), x(v, 2), x(v, 3)])
    for v in range(n):
        for c in range(4):
            for c2 in range(c + 1, 4):
                clauses.append([-x(v, c), -x(v, c2)])
    for (u, v) in sorted(E):
        for c in range(4):
            clauses.append([-x(u, c), -x(v, c)])
    fv = d["fixed_vertex"]
    clauses.append([x(fv, 0)])
    ncyc_cl = 0
    if not no_cycles:
        for cyc in d["cycles"]:
            m = len(cyc)
            assert m > 0 and len(set(cyc)) == m
            if m % 4 != 0:
                continue
            for s in range(4):
                clauses.append([-x(cyc[t], (s + t) % 4) for t in range(m)])
                ncyc_cl += 1
    nvars = 4 * n
    with open(out, "w") as fh:
        fh.write("p cnf %d %d\n" % (nvars, len(clauses)))
        for cl in clauses:
            fh.write(" ".join(map(str, cl)) + " 0\n")
    print("wrote %s: %d vars, %d clauses (vertices %d, edges %d, cycles %d, cycle clauses %d, no_cycles=%s)"
          % (out, nvars, len(clauses), n, len(E), 0 if no_cycles else len(d["cycles"]), ncyc_cl, no_cycles))


if __name__ == "__main__":
    main()
