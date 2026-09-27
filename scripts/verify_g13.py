"""Recheck alpha(G_13) = 36 from the certificates, with one command (notes/g13.md).

The claim rests on four facts:
  - the 15 sets KNOWN_G13 of scripts/experiments/largest_sets_whole_circles.py are independent sets of 36 points
    of G_13, so alpha(G_13) >= 36;
  - part A: the four formulas E37_A_c (c = 6, 7, 9, 11) of scripts/g13/enum_cert.py are unsatisfiable
    (certificates/g13_alpha_part_a_checks.txt);
  - part B: formula E37_B of scripts/g13/enum_cert.py is unsatisfiable under every leaf of
    certificates/g13_alpha_part_b_cubes.icnf (certificates/g13_alpha_part_b_checks.txt.gz);
  - the leaves of that file cover every assignment.
By the argument of notes/g13.md, section 3, parts A and B give alpha(G_13) <= 36.

This script checks the stored 6-colouring of G_13 (data/small_plane_colourings.json) and the 36-point sets, and that
the leaves cover every assignment. With scripts/g13/g13_audit.py present, it audits the formulas, from their text
alone. It rebuilds every formula and compares its SHA-256 with the one in the logs. It solves each formula with
kissat, which writes a DRAT proof, and has drat-trim check the proof. With --cake-lpr, drat-trim also writes the proof
in binary LRAT format, and cake_lpr checks that one too: cake_lpr is a proof checker verified in the HOL4 theorem
prover and compiled by the verified CakeML compiler. The proofs are deleted after the checks.

usage: python3 scripts/verify_g13.py --kissat PATH --drat-trim PATH [--cake-lpr PATH] [--jobs N]
                                     [--part a|b|all] [--sample K] [--shard K/N] [--from-log] [--log FILE]
       python3 scripts/verify_g13.py --no-solve
Formulas and proofs go to the directory in HN_OUT (default /tmp/hn), and one line per formula to the log (by default
g13_verify.txt there). A run that was interrupted resumes: the formulas that already have a confirmed line in the
log are skipped. With --sample K, only K leaves of part B, chosen at random, are solved: a spot check, not a proof.
With --no-solve, no solver runs: only the colouring, the independent sets, the audit, the cover and the SHA-256 of
every formula are checked. With --shard K/N, only the formulas of rank K, K + N, K + 2N, ... are solved, and the
cover only in shard 1. The logs of all N shards, concatenated, confirm every formula; a run over them with
--from-log solves nothing, fails if a formula has no confirmed line, and otherwise ends with CONFIRMED. The last
line says what was confirmed; the exit status is 0 only if everything that was checked is confirmed.
"""
import argparse
import gzip
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
HN = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "g13"))
sys.path.insert(0, os.path.join(HERE, "experiments"))
import g13                                               # noqa: E402
import enum_cert                                         # noqa: E402
from cuber2 import check_cover as covers                 # noqa: E402
from largest_sets_whole_circles import KNOWN_G13         # noqa: E402
try:
    import g13_audit as au                               # noqa: E402
except ImportError:
    au = None

CERT = os.path.join(HN, "certificates")
LOG_A = os.path.join(CERT, "g13_alpha_part_a_checks.txt")
LOG_B = os.path.join(CERT, "g13_alpha_part_b_checks.txt.gz")
CUBES = os.path.join(CERT, "g13_alpha_part_b_cubes.icnf")
COLOURINGS = os.path.join(HN, "data", "small_plane_colourings.json")
CIRCLES = (6, 7, 9, 11)
LINE = re.compile(r"(\S+): sha256 ([0-9a-f]{64}); kissat UNSAT in [\d.]+ s \(proof [\d.]+ MB\); "
                  r"drat-trim VERIFIED, \d+ of \d+ lemmas in core in [\d.]+ s$")


def formula_a(c):
    """E37_A_c: an independent dominating set of at least 37 points containing 0 and the whole circle N = c"""
    return enum_cert.build(37, 0, 0, [], False, c)[0].text()


def formula_b():
    """E37_B: an independent dominating set of at least 37 points, with no point together with its whole circle
    N = 6, 7, 9 or 11, lex-leader on the first 25 positions of lex_order() under the automorphisms"""
    return enum_cert.build(37, 0, 25, [], True, 0)[0].text()


def logged(path, opener=open):
    """name -> SHA-256, from a log of scripts/g13/certify.py (every line must say UNSAT and VERIFIED)"""
    out = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = LINE.match(line.rstrip("\n"))
            if not m or m.group(1) in out:
                raise ValueError(f"{path}: unexpected line {line!r}")
            out[m.group(1)] = m.group(2)
    return out


def read_cubes(path):
    """the leaves of the cube file, in order: 'a lit ... 0' (open after the cuber) or 'c closed lit ... 0' (closed
    by unit propagation); every one is a leaf of the tree, and every one is refuted"""
    out = []
    for line in open(path):
        if line.startswith("a "):
            out.append(tuple(map(int, line.split()[1:-1])))
        elif line.startswith("c closed "):
            out.append(tuple(map(int, line.split()[2:-1])))
    return out


def leaf_formula(base, cube):
    """formula E37_B (base: its header and the rest) with the literals of a leaf as unit clauses"""
    head, body = base
    nv, nc = head.split()[2:4]
    return f"p cnf {nv} {int(nc) + len(cube)}\n" + body + "".join(f"{l} 0\n" for l in cube)


class LeafDigests:
    """The SHA-256 of leaf_formula(base, cube) without writing the formula out: the hash of its header and of the
    clauses of E37_B is computed once for each number of unit clauses, then copied and extended by the cube."""

    def __init__(self, base):
        self.head, self.body = base
        self.prefix, self.lock = {}, threading.Lock()

    def __call__(self, cube):
        n = len(cube)
        with self.lock:
            if n not in self.prefix:
                nv, nc = self.head.split()[2:4]
                self.prefix[n] = hashlib.sha256(f"p cnf {nv} {int(nc) + n}\n{self.body}".encode())
            h = self.prefix[n].copy()
        h.update("".join(f"{l} 0\n" for l in cube).encode())
        return h.hexdigest()


def cover_proof(cubes):
    """The cover formula, whose clauses are the negated leaves, and a proof that it is unsatisfiable: every node of
    the decision tree, from the deepest up, gets the negation of its cube, resolved from its two children's
    clauses; the root's is the empty clause. Returns (DIMACS text, LRAT text, DRAT text)."""
    ids = {c: k + 1 for k, c in enumerate(cubes)}
    nv = max(abs(l) for c in cubes for l in c)
    cnf = f"p cnf {nv} {len(cubes)}\n" + "".join(" ".join(str(-l) for l in c) + " 0\n" for c in cubes)
    lrat, drat, top = [], [], [len(cubes)]

    def build(group, depth):
        if len(group) == 1 and len(group[0]) == depth:
            return ids[group[0]]
        if any(len(c) == depth for c in group) or len({abs(c[depth]) for c in group}) != 1:
            raise ValueError("the leaves are not the leaves of a binary decision tree")
        v = abs(group[0][depth])
        sides = [[c for c in group if c[depth] == v], [c for c in group if c[depth] == -v]]
        if not sides[0] or not sides[1]:
            raise ValueError("a node of the decision tree has only one child")
        a, b = build(sides[0], depth + 1), build(sides[1], depth + 1)
        top[0] += 1
        clause = " ".join(str(-l) for l in group[0][:depth])
        lrat.append(f"{top[0]} {clause} 0 {a} {b} 0\n".replace("  ", " "))
        drat.append(f"{clause} 0\n".lstrip())
        return top[0]

    build(list(cubes), 0)
    return cnf, "".join(lrat), "".join(drat)


def check_cover_proof(cubes, a):
    """check the cover formula's proof with cake_lpr (LRAT), or else with drat-trim (DRAT): (confirmed, report)"""
    cnf, lrat, drat = cover_proof(cubes)
    stem = os.path.join(a.out, "verify_g13_cover")
    try:
        with open(stem + ".cnf", "w") as fh:
            fh.write(cnf)
        with open(stem + ".proof", "w") as fh:
            fh.write(lrat if a.cake_lpr else drat)
        checker = a.cake_lpr if a.cake_lpr else a.drat_trim
        r = subprocess.run([checker, stem + ".cnf", stem + ".proof"], capture_output=True, text=True)
        ok = re.search(r"^s VERIFIED UNSAT" if a.cake_lpr else r"^s VERIFIED", r.stdout, re.M) is not None
        tool = "cake_lpr" if a.cake_lpr else "drat-trim"
        sha = hashlib.sha256(cnf.encode()).hexdigest()
        return ok, (f"the cover formula ({len(cubes)} negated leaves, sha256 {sha}) is unsatisfiable: "
                    f"{tool} {'VERIFIED' if ok else 'FAILED'}")
    finally:
        for ext in (".cnf", ".proof"):
            if os.path.exists(stem + ext):
                os.remove(stem + ext)


def solve(name, text, a):
    """kissat, drat-trim and, with --cake-lpr, cake_lpr on one formula: (confirmed, report)"""
    stem = os.path.join(a.out, f"verify_g13_{name}")
    with open(stem + ".cnf", "w") as fh:
        fh.write(text)
    try:
        t0 = time.time()
        k = subprocess.run([a.kissat, stem + ".cnf", stem + ".drat"], capture_output=True, text=True)
        t1 = time.time()
        report = f"kissat {'UNSAT' if k.returncode == 20 else f'exit status {k.returncode}'} in {t1 - t0:.0f} s"
        if k.returncode != 20:
            return False, report
        lrat = ["-L", stem + ".lrat", "-C"] if a.cake_lpr else []
        d = subprocess.run([a.drat_trim, stem + ".cnf", stem + ".drat", "-t", "200000"] + lrat,
                           capture_output=True, text=True)
        t2 = time.time()
        ok = re.search(r"^s VERIFIED", d.stdout, re.M) is not None
        report += f"; drat-trim {'VERIFIED' if ok else 'NOT VERIFIED'} in {t2 - t1:.0f} s"
        if not ok or not a.cake_lpr:
            return ok, report
        c = subprocess.run([a.cake_lpr, stem + ".cnf", stem + ".lrat"], capture_output=True, text=True)
        t3 = time.time()
        ok = re.search(r"^s VERIFIED UNSAT", c.stdout, re.M) is not None
        why = "" if ok else " (" + " ".join((c.stdout + c.stderr).split())[-120:] + ")"
        return ok, report + f"; cake_lpr {'VERIFIED UNSAT' if ok else 'FAILED'}{why} in {t3 - t2:.0f} s"
    finally:
        for ext in (".cnf", ".drat", ".lrat"):
            if os.path.exists(stem + ext):
                os.remove(stem + ext)


def main():
    ap = argparse.ArgumentParser(description="Recheck alpha(G_13) = 36 from the certificates.")
    ap.add_argument("--kissat")
    ap.add_argument("--drat-trim")
    ap.add_argument("--cake-lpr", default=None, help="also check each proof with cake_lpr (LRAT)")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--part", choices=("a", "b", "all"), default="all")
    ap.add_argument("--sample", type=int, default=0, help="solve only this many random leaves of part B")
    ap.add_argument("--seed", type=int, default=None, help="seed for --sample (default: random)")
    ap.add_argument("--shard", default=None, help="K/N: solve only the formulas of rank K, K + N, ... (1 <= K <= N)")
    ap.add_argument("--from-log", action="store_true",
                    help="solve nothing: every formula must have a confirmed line in the log (e.g. the shards' logs)")
    ap.add_argument("--no-solve", action="store_true", help="check the sets, the audit, the cover and hashes only")
    ap.add_argument("--log", default=None, help="default: g13_verify.txt in HN_OUT")
    a = ap.parse_args()
    if not a.no_solve and not (a.kissat and a.drat_trim):
        ap.error("solving needs --kissat and --drat-trim (or give --no-solve)")
    shard = None
    if a.shard:
        m = re.fullmatch(r"(\d+)/(\d+)", a.shard)
        if not m or not 1 <= int(m.group(1)) <= int(m.group(2)):
            ap.error("--shard takes K/N with 1 <= K <= N")
        shard = int(m.group(1)), int(m.group(2))
    a.out = os.environ.get("HN_OUT", "/tmp/hn")
    os.makedirs(a.out, exist_ok=True)
    log = a.log or os.path.join(a.out, "g13_verify.txt")
    failures, lock = [], threading.Lock()

    def say(line, ok=True):
        with lock:
            print(line, flush=True)
            if not ok:
                failures.append(line)
            if not a.no_solve:
                with open(log, "a") as fh:
                    fh.write(line + "\n")

    done = {}      # name -> SHA-256, for the formulas the log confirms (with cake_lpr, when it is asked for)
    if not a.no_solve and os.path.exists(log):
        for line in open(log):
            m = re.match(r"(E37_A\d+|leaf \d+): sha256 ([0-9a-f]{64}) as in the certificate; .*; confirmed$", line)
            if m and (not a.cake_lpr or "; cake_lpr VERIFIED UNSAT in " in line):
                done[m.group(1)] = m.group(2)
    if done:
        print(f"{log}: the formulas it confirms already ({len(done)}) are not solved again", flush=True)

    # the stored 6-colouring and the 36-point sets
    col = json.load(open(COLOURINGS))["G_q"]["13"]
    bad = [(u, v) for u, v in g13.EDGES if col[u] == col[v]]
    say(f"the stored colouring: {len(set(col))} colours, {len(bad)} of the {len(g13.EDGES)} edges of G_13 "
        f"monochromatic", ok=not bad and len(col) == 169 and set(col) == set(range(6)))
    adj = [set(n) for n in g13.ADJ]
    indep = [len(set(S)) == 36 and all(v not in adj[u] for u in S for v in S) for S in KNOWN_G13]
    say(f"the {len(KNOWN_G13)} known sets: 36 points each, independent: {all(indep)}",
        ok=all(indep) and len(KNOWN_G13) == 15)

    def audit(name, run_audit):
        if au is None:
            print(f"{name}: scripts/g13/g13_audit.py is missing, no audit", flush=True)
            return
        try:
            found = run_audit()
            say(f"{name}: the audit (scripts/g13/g13_audit.py) finds every clause as intended: {found}")
        except au.AuditError as e:
            say(f"{name}: the audit fails: {e}", ok=False)

    jobs = []
    if a.part in ("a", "all"):
        sha_a = logged(LOG_A)
        if sorted(sha_a) != sorted(f"E37_A{c}.cnf" for c in CIRCLES):
            say(f"part A: the log has the formulas {sorted(sha_a)}, not E37_A6, E37_A7, E37_A9, E37_A11", ok=False)
        for c in CIRCLES:
            audit(f"formula E37_A{c}", lambda c=c: au.audit_formula_a(c, formula_a(c)))
            jobs.append((f"E37_A{c}", lambda c=c: formula_a(c), sha_a.get(f"E37_A{c}.cnf"), None))
    if a.part in ("b", "all"):
        base_text = formula_b()
        audit("formula E37_B", lambda: au.audit_formula_b(base_text))
        cubes = read_cubes(CUBES)
        cover = covers(cubes)
        closed = sum(1 for line in open(CUBES) if line.startswith("c closed "))
        say(f"part B: {len(cubes)} leaves ({closed} closed by unit propagation); they cover every assignment: {cover}",
            ok=cover)
        if not a.no_solve and (shard is None or shard[0] == 1):
            ok, report = check_cover_proof(cubes, a)
            say(f"part B: {report}", ok=ok)
        sha_b = logged(LOG_B, gzip.open)
        if sorted(sha_b) != sorted(f"E37_B_leaf{i}" for i in range(len(cubes))):
            say(f"part B: the log does not have exactly one line per leaf ({len(sha_b)} lines)", ok=False)
        todo = list(range(len(cubes)))
        if a.sample:
            seed = a.seed if a.seed is not None else random.randrange(10 ** 9)
            todo = sorted(random.Random(seed).sample(todo, min(a.sample, len(todo))))
            print(f"part B: a random sample of {len(todo)} leaves (seed {seed})", flush=True)
        base = base_text.split("\n", 1)
        digest = LeafDigests(base)
        for i in todo:
            jobs.append((f"leaf {i}", lambda i=i: leaf_formula(base, cubes[i]), sha_b.get(f"E37_B_leaf{i}"),
                         lambda i=i: digest(cubes[i])))

    if shard:
        jobs = jobs[shard[0] - 1::shard[1]]
        print(f"shard {shard[0]} of {shard[1]}: {len(jobs)} formulas", flush=True)

    def run(job):
        name, build, sha, digest = job
        if a.no_solve or a.from_log or done.get(name) == sha:        # nothing to solve: hash without writing out
            got = digest() if digest else hashlib.sha256(build().encode()).hexdigest()
            if sha is None or got != sha:
                say(f"{name}: the formula does not match the certificate log", ok=False)
            elif not a.no_solve and done.get(name) != sha:
                say(f"{name}: no confirmed line in the log", ok=False)
            return
        text = build()
        if sha is None or hashlib.sha256(text.encode()).hexdigest() != sha:
            say(f"{name}: the formula does not match the certificate log", ok=False)
            return
        ok, report = solve(name.replace(" ", ""), text, a)
        say(f"{name}: sha256 {sha} as in the certificate; {report}; " + ("confirmed" if ok else "NOT CONFIRMED"),
            ok=ok)

    with ThreadPoolExecutor(max(1, a.jobs)) as ex:
        list(ex.map(run, jobs))

    what = "the colouring, the 36-point sets" + (", the audit of the formulas" if au else "")
    what += (", the cover" if a.part != "a" else "") + f" and the SHA-256 of {len(jobs)} formula"
    what += "s" if len(jobs) != 1 else ""
    if not a.no_solve:
        what += (", and every one of them unsatisfiable with a proof checked by drat-trim"
                 + (" and by cake_lpr" if a.cake_lpr else "") + (", as the log records" if a.from_log else ""))
    if failures:
        print(f"NOT CONFIRMED: {len(failures)} failures, the first: {failures[0]}")
        sys.exit(1)
    if a.no_solve or a.sample or a.part != "all" or shard or au is None:
        print(f"PARTIAL CHECK PASSED: {what}. This alone does not prove alpha(G_13) = 36.")
    else:
        print(f"CONFIRMED: {what}. So alpha(G_13) = 36.")


if __name__ == "__main__":
    main()
