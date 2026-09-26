"""Independent sets of the anisotropic plane G_17, towards alpha(G_17) <= 57 (notes/local_colourings.md §4).

G_17 = Cay(F_17^2, {N = 1}), N(x, y) = x^2 - 3y^2 (3 is the least non-square mod 17); the vertex (x, y) is
v = 17x + y. Every vertex has 18 neighbours, and the graph is vertex-transitive: translations and the maps
z -> l z + a, z -> l conj(z) + a (N(l) = 1) are automorphisms, 10 404 of them. Since 5 * 57 < 289, an
independence number of at most 57 would give chi(G_17) >= 6.

**Circles.** For c != 0 the circle C_c = {N = c} has 18 points. It contains no unit distance exactly for
c in {4, 5, 9, 11, 12, 14, 15}, the *independent circles*.

**The rosette.** kissat finds independent sets of 57 points, twice in the same orbit of the automorphism
group; ROSETTE below is one of them, translated so that 0 is its centre. It contains 0 and the whole circle
C_12, and 38 points on the circles N = 3, 4, 5, 6, 9, 14.

**Part A.** Let S be independent, let u be in S, and suppose the whole circle u + C_c lies in S. Then C_c is
an independent circle, and after translating u to 0 the other points of S lie in the region R_c of the
vertices adjacent to none of the 19 points of {0} + C_c. The formula A_c says that at least 58 - 19 = 39
vertices of R_c form an independent set. Every A_c is unsatisfiable (kissat, with a DRAT proof checked by
drat-trim; certificates/g17_part_a_checks.txt), so no independent set of 58 points contains a point
together with its whole circle.

**Part B, not settled.** The remaining case: an independent set of 58 points containing no point together
with its whole independent circle. formula_b() writes it, with vertex 0 in the set (translations preserve
the condition). If it is unsatisfiable, alpha(G_17) = 57 and chi(G_17) = 6.

The cardinality constraints use a totalizer: each node p of a balanced binary tree over the literals has
outputs o_j(p), j = 1..min(t, size), read 'at least j of the literals below p are true'. For children L, R
and a in [0, |L|], b in [0, |R|] with a + b + 1 <= min(t, |p|) there is the clause
    o_{a+1}(L) | o_{b+1}(R) | -o_{a+b+1}(p)     (the literal o_{|L|+1}(L), or o_{|R|+1}(R), is left out),
and the unit clause o_t(root). The assignment o_j(p) = [at least j true below p] satisfies every clause,
so the formula is satisfiable whenever at least t literals can be true (tests/test_g17.py checks both
directions exhaustively on small cases).

usage: python3 scripts/g17_alpha.py --kissat PATH --drat-trim PATH [--log FILE]      (part A)
       python3 scripts/g17_alpha.py --write-b FILE                                  (the formula of part B)
The formulas and proofs of part A are written to the directory in HN_OUT (default /tmp/hn), and the log
(by default g17_part_a_checks.txt there) gets one line per circle.
"""
import argparse
import hashlib
import os
import re
import subprocess
import time

Q, N0 = 17, 3
T = 58


def norm(z):
    return (z[0] * z[0] - N0 * z[1] * z[1]) % Q


def add(a, b):
    return ((a[0] + b[0]) % Q, (a[1] + b[1]) % Q)


def sub(a, b):
    return ((a[0] - b[0]) % Q, (a[1] - b[1]) % Q)


POINTS = [(x, y) for x in range(Q) for y in range(Q)]
INDEX = {z: Q * z[0] + z[1] for z in POINTS}
UNITS = [z for z in POINTS if norm(z) == 1]
EDGES = sorted({tuple(sorted((INDEX[z], INDEX[add(z, u)]))) for z in POINTS for u in UNITS})


def circle(c):
    return [z for z in POINTS if norm(z) == c]


def independent(points):
    return all(norm(sub(a, b)) != 1 for i, a in enumerate(points) for b in points[i + 1:])


INDEPENDENT_CIRCLES = [c for c in range(1, Q) if independent(circle(c))]

# the 57-point rosette, centred at 0: {0}, the circle N = 12, and 38 points of the circles N = 3, 4, 5, 6, 9, 14
ROSETTE = [
    (0, 0), (0, 4), (0, 8), (0, 9), (0, 13), (1, 6), (1, 11), (3, 0), (3, 4), (3, 8), (3, 9), (3, 13),
    (4, 2), (4, 15), (6, 4), (6, 5), (6, 8), (6, 9), (6, 12), (6, 13), (7, 1), (7, 2), (7, 3), (7, 6),
    (7, 7), (7, 10), (7, 11), (7, 14), (7, 15), (7, 16), (10, 1), (10, 2), (10, 3), (10, 6), (10, 7),
    (10, 10), (10, 11), (10, 14), (10, 15), (10, 16), (11, 5), (11, 12), (13, 2), (13, 3), (13, 7), (13, 10),
    (13, 14), (13, 15), (14, 0), (14, 4), (14, 8), (14, 9), (14, 13), (16, 2), (16, 6), (16, 11), (16, 15),
]


def totalizer(lits, t, top):
    """clauses saying that at least t of lits are true, with fresh variables from top + 1; (clauses, top)"""
    clauses = []

    def build(lo, hi):
        nonlocal top
        if hi - lo == 1:
            return [lits[lo]]
        mid = (lo + hi) // 2
        left, right = build(lo, mid), build(mid, hi)
        m = min(t, hi - lo)
        out = list(range(top + 1, top + m + 1))
        top += m
        for a in range(len(left) + 1):
            for b in range(len(right) + 1):
                j = a + b + 1
                if j > m:
                    continue
                cl = ([left[a]] if a < len(left) else []) + ([right[b]] if b < len(right) else [])
                clauses.append(cl + [-out[j - 1]])
        return out

    root = build(0, len(lits))
    clauses.append([root[t - 1]])
    return clauses, top


def dimacs(nvars, clauses):
    return f"p cnf {nvars} {len(clauses)}\n" + "".join(" ".join(map(str, c)) + " 0\n" for c in clauses)


def region(c):
    """the vertices adjacent to none of {0} + C_c, in the order of POINTS"""
    core = [(0, 0)] + circle(c)
    banned = set(core) | {add(p, u) for p in core for u in UNITS}
    return [z for z in POINTS if z not in banned]


def formula_a(c):
    """A_c: at least T - 19 vertices of the region R_c form an independent set; (text, |R_c|)"""
    reg = region(c)
    pos = {z: i + 1 for i, z in enumerate(reg)}
    clauses = [[-pos[z], -pos[add(z, u)]] for z in reg for u in UNITS if add(z, u) in pos and pos[add(z, u)] > pos[z]]
    card, top = totalizer(list(range(1, len(reg) + 1)), T - 19, len(reg))
    return dimacs(top, clauses + card), len(reg)


def formula_b():
    """B: an independent set of T points containing vertex 0 and no point together with its whole
    independent circle; s_v = 1 + v"""
    s = lambda z: 1 + INDEX[z]
    clauses = [[-(1 + u), -(1 + v)] for u, v in EDGES]
    for c in INDEPENDENT_CIRCLES:
        cc = circle(c)
        clauses += [[-s(u)] + [-s(add(u, w)) for w in cc] for u in POINTS]
    card, top = totalizer([s(z) for z in POINTS], T, len(POINTS))
    return dimacs(top, clauses + card + [[s((0, 0))]])


def check_a(c, kissat, drat_trim, out):
    text, size = formula_a(c)
    stem = os.path.join(out, f"g17_a{c}")
    with open(stem + ".cnf", "w") as fh:
        fh.write(text)
    sha = hashlib.sha256(text.encode()).hexdigest()
    t0 = time.time()
    r = subprocess.run([kissat, stem + ".cnf", stem + ".drat"], capture_output=True, text=True)
    t1 = time.time()
    verdict = {10: "SAT", 20: "UNSAT"}.get(r.returncode, "UNKNOWN")
    checked, core = "not run", ""
    if verdict == "UNSAT":
        d = subprocess.run([drat_trim, stem + ".cnf", stem + ".drat", "-t", "100000"], capture_output=True, text=True)
        text_out = d.stdout + d.stderr
        checked = "VERIFIED" if re.search(r"^s VERIFIED", text_out, re.M) else "NOT VERIFIED"
        m = re.search(r"(\d+) of (\d+) lemmas in core", text_out)
        core = f", {m.group(1)} of {m.group(2)} lemmas in core" if m else ""
    t2 = time.time()
    if os.path.exists(stem + ".drat"):
        os.remove(stem + ".drat")
    head = text.split("\n", 1)[0]
    return (f"A_{c}: region of {size} vertices, at least {T - 19} of them; {head}, sha256 {sha}; "
            f"kissat {verdict} in {t1 - t0:.0f} s; drat-trim {checked}{core} in {t2 - t1:.0f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kissat")
    ap.add_argument("--drat-trim")
    ap.add_argument("--log", default=None, help="default: g17_part_a_checks.txt in HN_OUT")
    ap.add_argument("--write-b", default=None, help="write the formula of part B to this file and stop")
    a = ap.parse_args()
    if a.write_b:
        text = formula_b()
        with open(a.write_b, "w") as fh:
            fh.write(text)
        print(f"{text.split(chr(10), 1)[0]}; sha256 {hashlib.sha256(text.encode()).hexdigest()}")
        return
    if not (a.kissat and a.drat_trim):
        ap.error("part A needs --kissat and --drat-trim")
    out = os.environ.get("HN_OUT", "/tmp/hn")
    os.makedirs(out, exist_ok=True)
    log = a.log or os.path.join(out, "g17_part_a_checks.txt")
    for c in INDEPENDENT_CIRCLES:
        line = check_a(c, a.kissat, a.drat_trim, out)
        print(line, flush=True)
        with open(log, "a") as fh:
            fh.write(line + "\n")


if __name__ == "__main__":
    main()
