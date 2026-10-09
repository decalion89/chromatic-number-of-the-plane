#!/usr/bin/env python3
"""Compare the referee CNF with the stored CNF as multisets of clauses,
after renaming the referee variables x(v,c) = c*n + v + 1 to the
vertex-major numbering 4*v + c + 1 (the layout the stored CNF appears to
use; the comparison itself decides whether that guess is right).

Usage: python3 compare_cnf.py REF.cnf STORED.cnf N_VERTICES
"""
import sys
from collections import Counter


def read_cnf(path):
    hdr = None
    clauses = []
    cur = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("c"):
                continue
            if s.startswith("p"):
                hdr = s.split()
                continue
            for tok in s.split():
                lit = int(tok)
                if lit == 0:
                    clauses.append(cur)
                    cur = []
                else:
                    cur.append(lit)
    assert not cur, "unterminated clause"
    return hdr, clauses


def main():
    ref, stored, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
    h1, c1 = read_cnf(ref)
    h2, c2 = read_cnf(stored)
    print("ref header", h1, "parsed clauses", len(c1))
    print("stored header", h2, "parsed clauses", len(c2))
    assert int(h1[3]) == len(c1) and int(h2[3]) == len(c2)

    def ren(lit):
        v = abs(lit) - 1
        c, vert = divmod(v, n)
        new = 4 * vert + c + 1
        return new if lit > 0 else -new

    A = Counter(tuple(sorted(set(ren(l) for l in cl))) for cl in c1)
    B = Counter(tuple(sorted(set(cl))) for cl in c2)
    print("distinct clauses: ref %d, stored %d" % (len(A), len(B)))
    onlyA = A - B
    onlyB = B - A
    print("clauses only in ref (after renaming):", sum(onlyA.values()))
    print("clauses only in stored:", sum(onlyB.values()))
    for cl in list(onlyA)[:5]:
        print("  ref-only example:", cl)
    for cl in list(onlyB)[:5]:
        print("  stored-only example:", cl)
    # also: any tautologies / duplicate literals in stored?
    taut = sum(1 for cl in c2 if any(-l in cl for l in cl))
    dupl = sum(1 for cl in c2 if len(set(cl)) != len(cl))
    maxvar = max(abs(l) for cl in c2 for l in cl)
    print("stored: tautological clauses %d, clauses with repeated literal %d, max var %d" % (taut, dupl, maxvar))
    same = (A == B)
    print("IDENTICAL AS CLAUSE MULTISETS UNDER RENAMING:", same)
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
