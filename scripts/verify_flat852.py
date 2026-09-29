"""Recheck that data/flat852/core852a.json is a 5-chromatic unit-distance graph (notes/flat852.md).

The claim rests on four facts:
  - every vertex is a point of Q(zeta21), given by 12 integers c_0..c_11 as (c_0 + c_1 z + ... + c_11 z^11)/7 with
    z = zeta21, and the vertices are distinct;
  - every edge joins two vertices whose difference d is one of the 126 unit vectors U = mu42 u omega*mu42 u
    conj(omega)*mu42 listed in the file, and d*conj(d) = 1 holds exactly in Q(zeta21);
  - the stored 5-colouring is proper, so the graph is 5-colourable;
  - the formula "the graph is 4-colourable, with the colours of one triangle fixed" is unsatisfiable.
The first three are checked here in exact integer arithmetic, with no library. For the last one, the script checks
that the stored CNF (data/flat852/core852a.cnf, the formula of the certificate in data/flat852/logs/) is exactly that
formula for this graph; with --kissat and --drat-trim it also writes the formula again in its own encoding, solves it
with kissat, which writes a DRAT proof, and has drat-trim check the proof (about 20 and 35 minutes on one core).

usage: python3 scripts/verify_flat852.py [--graph FILE] [--colouring FILE] [--cnf FILE]
                                         [--kissat PATH --drat-trim PATH] [--workdir DIR]
The last line begins with CONFIRMED, or with NOT CONFIRMED and the first failure; the exit status is 0 only if
everything that was checked is confirmed."""
import argparse
import itertools
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N = 21                                   # z = zeta21; elements of Q(zeta21) are 12 integers over 7


def poly_divmod(a, b):
    """quotient and remainder of integer polynomials (lists, lowest degree first), b monic"""
    a = list(a)
    q = [0] * max(len(a) - len(b) + 1, 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1]
        q[i] = c
        for j, bj in enumerate(b):
            a[i + j] -= c * bj
    return q, a[:len(b) - 1]


def cyclotomic21():
    """Phi_21 = (x^21 - 1)(x - 1) / ((x^7 - 1)(x^3 - 1)), degree 12"""
    def mul(p, r):
        out = [0] * (len(p) + len(r) - 1)
        for i, pi in enumerate(p):
            for j, rj in enumerate(r):
                out[i + j] += pi * rj
        return out
    num = mul([-1] + [0] * 20 + [1], [-1, 1])
    den = mul([-1] + [0] * 6 + [1], [-1, 0, 0, 1])
    q, r = poly_divmod(num, den)
    assert not any(r) and len(q) == 13 and q[-1] == 1
    return q


PHI = cyclotomic21()


def reduce(p):
    return poly_divmod(list(p) + [0] * max(0, 13 - len(p)), PHI)[1] if len(p) > 12 else list(p) + [0] * (12 - len(p))


def mul(a, b):
    out = [0] * 23
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return reduce(out)


def conj(a):
    """z -> z^20 = z^-1"""
    out = [0] * 21
    for k, ak in enumerate(a):
        out[(N - k) % N] += ak
    return reduce(out)


def is_unit_vector(d):
    """(d/7) * conj(d/7) = 1, i.e. d * conj(d) = 49"""
    return mul(d, conj(d)) == [49] + [0] * 11


def load(path):
    with open(path) as fh:
        return json.load(fh)


def check_directions(U):
    """the 126 unit vectors: distinct, each of norm 1, closed under negation and conjugation"""
    S = {tuple(u) for u in U}
    return (len(U) == 126 and len(S) == 126 and all(is_unit_vector(u) for u in U)
            and all(tuple(-x for x in u) in S for u in U) and all(tuple(conj(u)) in S for u in U))


def check_graph(g):
    """(ok, message): distinct vertices, every edge a unit vector of U, the fixed triangle a triangle"""
    V = [tuple(v["exact7"]) for v in g["vertices"]]
    if len(set(V)) != len(V):
        return False, "two vertices coincide"
    U = [tuple(u) for u in g["U"]]
    if not check_directions(U):
        return False, "the list U is not 126 distinct unit vectors closed under negation and conjugation"
    for a, b, j in g["edges"]:
        d = tuple(y - x for x, y in zip(V[a], V[b]))
        if d != U[j]:
            return False, f"edge ({a}, {b}): its vector is not U[{j}]"
    if len({frozenset((a, b)) for a, b, _ in g["edges"]}) != len(g["edges"]):
        return False, "an edge is listed twice"
    adj = {frozenset((a, b)) for a, b, _ in g["edges"]}
    t = g["triangle_fixed"]
    if not all(frozenset(p) in adj for p in itertools.combinations(t, 2)):
        return False, f"the fixed vertices {t} are not a triangle"
    return True, f"{len(V)} distinct vertices, {len(g['edges'])} edges, every edge one of the 126 unit vectors"


def check_colouring(g, colours):
    n = len(g["vertices"])
    if len(colours) != n:
        return False, "the colouring has the wrong length"
    if any(colours[a] == colours[b] for a, b, _ in g["edges"]):
        return False, "the colouring is not proper"
    return True, f"a proper colouring with {len(set(colours))} colours"


def cnf_clauses(g, var):
    """the clauses of 'the graph is 4-colourable, the triangle coloured 0, 1, 2' for a variable map var(v, c)"""
    n = len(g["vertices"])
    cl = [tuple(var(v, c) for c in range(4)) for v in range(n)]
    cl += [(-var(a, c), -var(b, c)) for a, b, _ in g["edges"] for c in range(4)]
    cl += [(var(v, i),) for i, v in enumerate(g["triangle_fixed"])]
    return cl


def read_cnf(path):
    nv, clauses, cur = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("c"):
                continue
            if line.startswith("p"):
                nv = int(line.split()[2])
                continue
            for tok in line.split():
                x = int(tok)
                if x == 0:
                    clauses.append(tuple(cur))
                    cur = []
                else:
                    cur.append(x)
    return nv, clauses


def check_stored_cnf(g, path):
    """the stored CNF is exactly the formula of this graph in the encoding 4v + c + 1 (as a multiset of clauses)"""
    nv, clauses = read_cnf(path)
    want = cnf_clauses(g, lambda v, c: 4 * v + c + 1)
    norm = lambda cls: sorted(tuple(sorted(c)) for c in cls)
    if nv != 4 * len(g["vertices"]) or norm(clauses) != norm(want):
        return False, f"{os.path.relpath(path, ROOT)} is not the 4-colouring formula of this graph"
    return True, f"{os.path.relpath(path, ROOT)} is the 4-colouring formula of this graph ({len(clauses)} clauses)"


def solve(g, kissat, drat_trim, workdir):
    """own encoding (c*n + v + 1), kissat with a DRAT proof, drat-trim"""
    n = len(g["vertices"])
    cl = cnf_clauses(g, lambda v, c: c * n + v + 1)
    os.makedirs(workdir, exist_ok=True)
    cnf, proof = os.path.join(workdir, "flat852_4col.cnf"), os.path.join(workdir, "flat852_4col.drat")
    with open(cnf, "w") as fh:
        fh.write(f"p cnf {4 * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    k = subprocess.run([kissat, cnf, proof], capture_output=True, text=True)
    if "s UNSATISFIABLE" not in k.stdout:
        return False, "kissat did not answer UNSATISFIABLE"
    d = subprocess.run([drat_trim, cnf, proof, "-t", "20000"], capture_output=True, text=True)
    os.remove(proof)
    if "s VERIFIED" not in d.stdout:
        return False, "drat-trim did not verify the proof"
    return True, "kissat UNSATISFIABLE and drat-trim VERIFIED on a formula written again here"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", default=os.path.join(ROOT, "data/flat852/core852a.json"))
    ap.add_argument("--colouring", default=os.path.join(ROOT, "data/flat852/core852.5colouring.json"))
    ap.add_argument("--cnf", default=os.path.join(ROOT, "data/flat852/core852a.cnf"))
    ap.add_argument("--kissat")
    ap.add_argument("--drat-trim")
    ap.add_argument("--workdir", default=os.path.join(tempfile.gettempdir(), "flat852"))
    a = ap.parse_args()
    g = load(a.graph)
    col = load(a.colouring)
    col = col.get("colours", col.get("colouring")) if isinstance(col, dict) else col
    steps = [lambda: check_graph(g), lambda: check_colouring(g, col), lambda: check_stored_cnf(g, a.cnf)]
    if a.kissat and a.drat_trim:
        steps.append(lambda: solve(g, a.kissat, a.drat_trim, a.workdir))
    for step in steps:
        ok, msg = step()
        print(("ok: " if ok else "FAILED: ") + msg)
        if not ok:
            print("NOT CONFIRMED: " + msg)
            return 1
    solved = " and not 4-colourable (kissat, drat-trim)" if a.kissat and a.drat_trim else \
        "; not 4-colourable by the stored certificate (run with --kissat and --drat-trim to solve again)"
    print(f"CONFIRMED: {os.path.relpath(a.graph, ROOT)} is a unit-distance graph with a proper 5-colouring{solved}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
