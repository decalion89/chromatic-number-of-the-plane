"""Recheck the four-chromatic unit-distance graphs over real quadratic fields of notes/quadratic_planes.md.

For each field Q(sqrt d) in data/quadratic_planes/ the file qD.json holds a graph G whose vertices are points
((a + b sqrt d)/D, (c + e sqrt d)/D) of the plane over Q(sqrt d), given by the integers [a, b, c, e]. The checks:
  - the points are distinct, and every edge has length exactly 1: with da, db, dc, de the differences,
    (da + db sqrt d)^2 + (dc + de sqrt d)^2 = D^2, i.e. da^2 + d db^2 + dc^2 + d de^2 = D^2 and da db + dc de = 0;
  - the edge list is every pair of points at distance 1 (so G is the unit-distance graph on its points), and G has
    no triangle;
  - the stored 4-colouring is proper, so chi(G) <= 4;
  - for every vertex v the stored 3-colouring of G - v is proper, so G is vertex-critical if chi(G) = 4;
  - the stored CNF (qD.cnf) is exactly the formula "G is 3-colourable, with the edge fixed[0], fixed[1] coloured 0, 1"
    in the encoding 3v + c + 1. Its unsatisfiability is the certificate in qD.logs/ (kissat, drat-trim); with
    --kissat PATH --drat-trim PATH this script writes the formula again in its own encoding (c n + v + 1), solves it
    and checks the proof; with --cake-lpr PATH as well, drat-trim also writes the proof in LRAT form and cake_lpr,
    a proof checker verified in CakeML, checks it too, and the same is done for the stored formula qD.cnf itself
    (the hypothesis of the Lean files lean/Sqrt{d}.lean that do not check an LRAT proof). Then chi(G) = 4 and
    chi(Q(sqrt d)^2) >= 4.
  - the upper bound: when d is a nonzero square modulo p = 7 or p = 11 (both 3 mod 4), reduction modulo a prime of
    norm p maps the plane over Q(sqrt d) into the unit-distance graph of F_p^2 (Moorhouse, Lemma 8.2), and
    finite_planes.json holds a proper 4-colouring of F_7^2 and a proper 5-colouring of F_11^2, checked here.
usage: python3 scripts/verify_quadratic_planes.py [--field d ...] [--kissat PATH --drat-trim PATH [--cake-lpr PATH]]
       [--workdir DIR]
The last line begins with CONFIRMED or with NOT CONFIRMED; the exit status is 0 only if everything checked holds."""
import argparse
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "quadratic_planes")


def load(path):
    with open(path) as fh:
        return json.load(fh)


def unit(d, D, p, q):
    da, db, dc, de = (y - x for x, y in zip(p, q))
    return da * da + d * db * db + dc * dc + d * de * de == D * D and da * db + dc * de == 0


def check_graph(g):
    d, D = g["d"], g["D"]
    P = [tuple(p) for p in g["points"]]
    n = len(P)
    if len(set(P)) != n:
        return False, "two points coincide"
    E = {(min(a, b), max(a, b)) for a, b in g["edges"]}
    if len(E) != len(g["edges"]):
        return False, "an edge is listed twice"
    for a, b in E:
        if not unit(d, D, P[a], P[b]):
            return False, f"edge ({a}, {b}) does not have length 1"
    allu = {(a, b) for a in range(n) for b in range(a + 1, n) if unit(d, D, P[a], P[b])}
    if allu != E:
        return False, f"{len(allu - E)} pairs at distance 1 are not listed as edges"
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b); adj[b].add(a)
    if any(adj[a] & adj[b] for a, b in E):
        return False, "the graph has a triangle"
    return True, f"{n} distinct points of the plane over Q(sqrt{d}), {len(E)} edges, all pairs at distance exactly 1, no triangle"


def digits(s):
    return [-1 if ch == "-" else int(ch) for ch in s] if isinstance(s, str) else s


def check_colouring(g, col, k, drop=None):
    col = digits(col)
    n = len(g["points"])
    if drop is not None and (len(col) != n or col[drop] != -1):
        return False
    if len(col) != n or any(not (0 <= col[v] < k) for v in range(n) if v != drop):
        return False
    return all(col[a] != col[b] for a, b in g["edges"] if drop not in (a, b))


def check_critical(g):
    n = len(g["points"])
    crit = g["critical_3_colourings"]
    if sorted(crit, key=int) != [str(v) for v in range(n)]:
        return False, "the vertex-critical certificate does not have one colouring per vertex"
    for v in range(n):
        if not check_colouring(g, crit[str(v)], 3, drop=v):
            return False, f"the colouring for vertex {v} is not a proper 3-colouring of G minus {v}"
    return True, f"for each of the {n} vertices v, a proper 3-colouring of G minus v"


def clauses(g, var):
    n = len(g["points"])
    cl = [tuple(var(v, c) for c in range(3)) for v in range(n)]
    cl += [(-var(a, c), -var(b, c)) for a, b in g["edges"] for c in range(3)]
    a0, b0 = g["fixed_edge"]
    return cl + [(var(a0, 0),), (var(b0, 1),)]


def read_cnf(path):
    nv, out, cur = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("c"):
                continue
            if line.startswith("p"):
                nv = int(line.split()[2]); continue
            for tok in line.split():
                x = int(tok)
                if x == 0:
                    out.append(tuple(cur)); cur = []
                else:
                    cur.append(x)
    return nv, out


def check_cnf(g, path):
    nv, cl = read_cnf(path)
    norm = lambda cs: sorted(tuple(sorted(c)) for c in cs)
    a0, b0 = g["fixed_edge"]
    if (min(a0, b0), max(a0, b0)) not in {(min(a, b), max(a, b)) for a, b in g["edges"]}:
        return False, "the fixed pair is not an edge"
    if nv != 3 * len(g["points"]) or norm(cl) != norm(clauses(g, lambda v, c: 3 * v + c + 1)):
        return False, f"{os.path.basename(path)} is not the 3-colouring formula of this graph"
    return True, f"{os.path.basename(path)} is the 3-colouring formula of this graph ({len(cl)} clauses)"


def refute(cnf, kissat, drat_trim, workdir, cake_lpr, where, name):
    """kissat on the formula in the file cnf, drat-trim on its proof, and, with cake_lpr, cake_lpr on the proof in
    LRAT form; where names the formula in the message, and the proofs are workdir/name.drat and name.lrat."""
    os.makedirs(workdir, exist_ok=True)
    proof, lrat = os.path.join(workdir, name + ".drat"), os.path.join(workdir, name + ".lrat")
    try:
        k = subprocess.run([kissat, cnf, proof], capture_output=True, text=True)
        if "s UNSATISFIABLE" not in k.stdout:
            return False, "kissat did not answer UNSATISFIABLE"
        r = subprocess.run([drat_trim, cnf, proof, "-t", "20000"] + (["-L", lrat, "-C"] if cake_lpr else []),
                           capture_output=True, text=True)
        if "s VERIFIED" not in r.stdout:
            return False, "drat-trim did not verify the proof"
        if not cake_lpr:
            return True, f"kissat UNSATISFIABLE and drat-trim VERIFIED on {where}"
        c = subprocess.run([cake_lpr, cnf, lrat], capture_output=True, text=True)
        if "s VERIFIED UNSAT" not in c.stdout:
            return False, "cake_lpr did not verify the LRAT proof"
        return True, f"kissat UNSATISFIABLE, drat-trim VERIFIED and cake_lpr VERIFIED UNSAT (LRAT) on {where}"
    finally:
        for f in (proof, lrat):
            if os.path.exists(f):
                os.remove(f)


def solve(g, kissat, drat_trim, workdir, cake_lpr=None):
    n = len(g["points"])
    cl = clauses(g, lambda v, c: c * n + v + 1)
    os.makedirs(workdir, exist_ok=True)
    cnf = os.path.join(workdir, f"q{g['d']}.cnf")
    with open(cnf, "w") as fh:
        fh.write(f"p cnf {3 * n} {len(cl)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in cl))
    return refute(cnf, kissat, drat_trim, workdir, cake_lpr, "a formula written again here", f"q{g['d']}")


def plane_colouring_ok(p, col):
    S = [(a, b) for a in range(p) for b in range(p) if (a * a + b * b) % p == 1]
    if len(col) != p * p:
        return False
    for x in range(p):
        for y in range(p):
            for a, b in S:
                if col[x * p + y] == col[((x + a) % p) * p + (y + b) % p]:
                    return False
    return True


def residue(d, p):
    return d % p and pow(d % p, (p - 1) // 2, p) == 1


def reduces(d, p):
    """the hypothesis of Moorhouse's Lemma 8.2 for p = 3 mod 4: d = 0 or a nonzero square mod p"""
    return d % p == 0 or residue(d, p)


def upper_bound(d, planes):
    """the best of: reduction at 7 (4 colours), the place over 2 when d = 3 mod 8 (4), reduction at 11 (5)"""
    def by(p):
        pl = planes[str(p)]
        if not plane_colouring_ok(p, pl["colouring"]):
            return None, f"the stored colouring of F_{p}^2 is not proper"
        how = "0" if d % p == 0 else "a nonzero square"
        return pl["colours"], f"d = {d} is {how} mod {p} and p = {p} is 3 mod 4: chi <= chi(F_{p}^2) <= " \
                              f"{pl['colours']} (Moorhouse, Lemma 8.2; stored colouring checked)"
    if reduces(d, 7):
        return by(7)
    if d % 8 == 3:
        return 4, f"d = {d} is 3 mod 8: chi <= 4 (Fischer 1990, Thm 10; notes/local_colourings.md, Proposition A)"
    if reduces(d, 11):
        return by(11)
    return None, "no upper bound below 7 recorded"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", type=int, nargs="*")
    ap.add_argument("--kissat")
    ap.add_argument("--drat-trim")
    ap.add_argument("--cake-lpr", help="also check each proof, in LRAT form, with cake_lpr")
    ap.add_argument("--workdir", default=os.path.join(tempfile.gettempdir(), "quadratic_planes"))
    a = ap.parse_args()
    planes = load(os.path.join(DATA, "finite_planes.json"))["planes"]
    ds = a.field or sorted(int(f[1:-5]) for f in os.listdir(DATA) if f.startswith("q") and f.endswith(".json"))
    summary = []
    for d in ds:
        g = load(os.path.join(DATA, f"q{d}.json"))
        steps = [lambda: check_graph(g),
                 lambda: (check_colouring(g, g["four_colouring"], 4), "a proper 4-colouring"),
                 lambda: check_critical(g),
                 lambda: check_cnf(g, os.path.join(DATA, f"q{d}.cnf"))]
        if a.kissat and a.drat_trim:
            steps.append(lambda: solve(g, a.kissat, a.drat_trim, a.workdir, a.cake_lpr))
        if a.kissat and a.drat_trim and a.cake_lpr:
            # the stored formula itself, the hypothesis of the Lean files of the fields without lrat_proof
            steps.append(lambda: refute(os.path.join(DATA, f"q{d}.cnf"), a.kissat, a.drat_trim, a.workdir, a.cake_lpr,
                                        f"q{d}.cnf itself", f"q{d}_stored"))
        for step in steps:
            ok, msg = step()
            print(f"Q(sqrt{d}): " + ("ok: " if ok else "FAILED: ") + msg)
            if not ok:
                print(f"NOT CONFIRMED: Q(sqrt{d}): {msg}")
                return 1
        ub, msg = upper_bound(d, planes)
        print(f"Q(sqrt{d}): upper bound: {msg}")
        if ub is None and "not proper" in msg:
            print(f"NOT CONFIRMED: Q(sqrt{d}): {msg}")
            return 1
        summary.append(f"chi(Q(sqrt{d})^2) {'= 4' if ub == 4 else '>= 4' if ub is None else f'in [4, {ub}]'}")
    how = ("kissat, drat-trim and cake_lpr run here" if a.cake_lpr else "kissat and drat-trim run here") \
        if a.kissat and a.drat_trim else \
        "non-3-colourability by the stored certificates (run with --kissat and --drat-trim to solve again)"
    print("CONFIRMED: " + "; ".join(summary) + f" ({how})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
