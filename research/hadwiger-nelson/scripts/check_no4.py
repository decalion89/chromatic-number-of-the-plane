"""Recheck, from the exact coordinates, that graphs in data/ have no proper 4-colouring.

usage: check_no4.py --kissat PATH --drat-trim PATH [--time SECONDS] [--log FILE] [FILE ...]

With no FILE, every graph in NO_FOUR below is checked, smallest first: the graphs of data/ that
data/README.md describes as having no proper 4-colouring (or as 5-chromatic). For each one:

1. the unit-distance graph is rebuilt from the exact coordinates by hn.graph.build_graph, which
   keeps a pair only if its squared distance is exactly 1; the file must have no repeated point;
2. the 4-colouring formula is written as in certificates/README.md, step 3: x(v, c) = 1 + 4v + c,
   one clause per vertex, then one clause per edge and colour, edges in lexicographic order of
   (u, v) with u < v. Then come the 12 unit clauses that pin the lexicographically first triangle
   (a, b, c) to colours 0, 1, 2, as step 4 does for de Grey's graph. The pin is sound: in a proper
   colouring the three vertices of a triangle receive distinct colours, and a permutation of the
   colours makes them 0, 1 and 2;
3. kissat solves the formula and writes a binary DRAT proof, and drat-trim checks the proof
   against the formula. The proof is then deleted; the formula's SHA-256 is recorded, so that the
   check can be repeated exactly.

One line per graph is appended to the log. The log (by default data_no4_checks.txt), the formula and
the proof are written to the directory in HN_OUT (default /tmp/hn); certificates/data_no4_checks.txt
is the record of the run of 26 September 2026.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from fractions import Fraction as Fr

HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HN_DIR)
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

NO_FOUR = [
    "five_247_c.json", "hunt_811.json", "five_247_b.json", "disc_1099.json", "disc_1113.json",
    "shrunk_1135.json", "five_247.json", "five_247_c_lambda.json", "five_tuned_1_1.json",
    "five_rho7.json", "five_tuned_1_3.json", "five_tuned_16_1_3_7_11.json", "five_tuned_2_1.json",
    "five_tuned_4_1.json", "five_twotune_small.json", "tight_hexagon_4159.json",
    "five_dense_2.json", "five_23.json", "five_symmetric.json", "five_23_v7.json",
    "five_order6_free.json", "five_twotune.json", "five_dense_10.json", "free_angle.json",
    "rarest_v139.json", "free_angle_both.json", "rarest_v199.json",
]
K = 4


def load(path):
    with open(path) as fh:
        d = json.load(fh)
    F = Field(tuple(d["field_generators"]))
    pts = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y]))
           for x, y in d["points"]]
    return d, pts


def first_triangle(n, adj):
    for a in range(n):
        for b in sorted(v for v in adj[a] if v > a):
            common = sorted(v for v in adj[a] & adj[b] if v > b)
            if common:
                return a, b, common[0]
    return None


def write_formula(n, edges, tri, path):
    def x(v, c):
        return 1 + K * v + c
    lines = [" ".join(str(x(v, c)) for c in range(K)) + " 0" for v in range(n)]
    for u, v in edges:
        lines.extend(f"{-x(u, c)} {-x(v, c)} 0" for c in range(K))
    if tri:
        for colour, v in enumerate(tri):
            lines.extend(f"{x(v, c) if c == colour else -x(v, c)} 0" for c in range(K))
    text = f"p cnf {n * K} {len(lines)}\n" + "\n".join(lines) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    return len(lines), hashlib.sha256(text.encode()).hexdigest()


def check(name, kissat, drat_trim, seconds, out):
    t0 = time.time()
    d, pts = load(os.path.join(HN_DIR, "data", name))
    g = build_graph(pts)
    assert g.n == len(pts), f"{name}: {len(pts) - g.n} repeated points"
    edges = sorted(g.edges())
    if "m" in d:
        assert len(edges) == d["m"], f"{name}: {len(edges)} edges, the file says {d['m']}"
    tri = first_triangle(g.n, g.adj)
    stem = os.path.join(out, name[:-5] + "_k4")
    ncl, sha = write_formula(g.n, edges, tri, stem + ".cnf")
    t1 = time.time()
    r = subprocess.run([kissat, f"--time={seconds}", stem + ".cnf", stem + ".drat"],
                       capture_output=True, text=True)
    t2 = time.time()
    verdict = {10: "SAT", 20: "UNSAT"}.get(r.returncode, "UNKNOWN")
    checked, core = "not run", ""
    if verdict == "UNSAT":
        c = subprocess.run([drat_trim, stem + ".cnf", stem + ".drat", "-t", "200000"],
                           capture_output=True, text=True)
        text = c.stdout + c.stderr
        checked = "VERIFIED" if re.search(r"^s VERIFIED", text, re.M) else "NOT VERIFIED"
        m = re.search(r"(\d+) of (\d+) lemmas in core", text)
        core = f"{m.group(1)} of {m.group(2)} lemmas in core" if m else ""
    t3 = time.time()
    if os.path.exists(stem + ".drat"):
        os.remove(stem + ".drat")
    return (f"{name}: {g.n} vertices, {len(edges)} edges; triangle {tri} pinned; "
            f"{ncl} clauses, sha256 {sha}; kissat {verdict} in {t2 - t1:.0f} s; "
            f"drat-trim {checked}{', ' + core if core else ''} in {t3 - t2:.0f} s; "
            f"graph built in {t1 - t0:.0f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--kissat", required=True)
    ap.add_argument("--drat-trim", required=True)
    ap.add_argument("--time", type=int, default=7200, help="kissat's limit per graph, in seconds")
    ap.add_argument("--log", default=None, help="default: data_no4_checks.txt in HN_OUT")
    a = ap.parse_args()
    out = os.environ.get("HN_OUT", "/tmp/hn")
    os.makedirs(out, exist_ok=True)
    log = a.log or os.path.join(out, "data_no4_checks.txt")
    for name in a.files or NO_FOUR:
        line = check(os.path.basename(name), a.kissat, a.drat_trim, a.time, out)
        print(line, flush=True)
        with open(log, "a") as fh:
            fh.write(line + "\n")


if __name__ == "__main__":
    main()
