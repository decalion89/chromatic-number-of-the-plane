"""Binary structured cubes: branch on the next unassigned variable of a fixed list (both signs), after unit
propagation.  A branch whose propagation fails is a closed leaf; a branch reaching MAXDEC decisions, or where
every listed variable is assigned, is an open leaf.  Every model of the formula satisfies exactly one leaf
(the tree branches on one variable both ways at every internal node); check_cover() verifies this from the
leaves alone.

usage: cuber2.py CNF MAXDEC OUT.icnf --vars v1,v2,...  [--solver m22] [--minpos P]
Output lines: 'a lit ... 0' (open leaf) or 'c closed lit ... 0'.
"""
import argparse
import os
import sys
from pysat.formula import CNF
from pysat.solvers import Solver


def build(clauses, maxdec, vars_, solver="m22", limit=10 ** 7, minpos=0):
    """minpos > 0: a node is also a leaf once minpos of its decisions are positive (for 'v in C0' variables:
    leaves then fix about the same number of points of C0, which balances their difficulty)"""
    s = Solver(name=solver, bootstrap_with=clauses)
    leaves = []

    def rec(path):
        if len(leaves) > limit:
            raise RuntimeError("too many leaves")
        ok, props = s.propagate(assumptions=path)
        if not ok:
            leaves.append(("closed", path))
            return
        if len(path) == maxdec or (minpos and sum(1 for l in path if l > 0) >= minpos):
            leaves.append(("open", path))
            return
        assigned = {abs(l) for l in props} | {abs(l) for l in path}
        nxt = next((v for v in vars_ if v not in assigned), None)
        if nxt is None:
            leaves.append(("open", path))
            return
        rec(path + [nxt])
        rec(path + [-nxt])

    rec([])
    s.delete()
    return leaves


def check_cover(paths, depth=0):
    """True when the paths are the leaves of a binary tree branching on one variable both ways at each node"""
    if len(paths) == 1 and len(paths[0]) == depth:
        return True
    if any(len(p) == depth for p in paths):
        return False
    lits = {p[depth] for p in paths}
    if len(lits) != 2 or len({abs(l) for l in lits}) != 1:
        return False
    v = abs(next(iter(lits)))
    return (check_cover([p for p in paths if p[depth] == v], depth + 1)
            and check_cover([p for p in paths if p[depth] == -v], depth + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cnf")
    ap.add_argument("maxdec", type=int)
    ap.add_argument("out")
    ap.add_argument("--vars", required=True)
    ap.add_argument("--solver", default="m22")
    ap.add_argument("--minpos", type=int, default=0, help="stop once this many decisions are positive")
    a = ap.parse_args()
    f = CNF(from_file=a.cnf)
    vars_ = [int(t) for t in a.vars.split(",")]
    leaves = build(f.clauses, a.maxdec, vars_, a.solver, minpos=a.minpos)
    sys.setrecursionlimit(100000)
    if not check_cover([tuple(p) for _, p in leaves]):
        sys.exit("cover check failed")
    with open(a.out + ".tmp", "w") as fh:          # written whole, then renamed: never a truncated cube file
        for kind, p in leaves:
            fh.write(("a " if kind == "open" else "c closed ") + " ".join(map(str, p)) + " 0\n")
    os.replace(a.out + ".tmp", a.out)
    n_open = sum(1 for k, _ in leaves if k == "open")
    print(f"{a.out}: {len(leaves)} leaves, {n_open} open, {len(leaves) - n_open} closed; cover OK")


if __name__ == "__main__":
    main()
