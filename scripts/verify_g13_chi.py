"""Recheck chi(G_13) = 6 from the certificates, with one command (notes/g13_chi.md).

The claim rests on four facts:
  - data/small_plane_colourings.json holds a proper 6-colouring of G_13, so chi(G_13) <= 6;
  - alpha(G_13) = 36 (notes/g13.md; scripts/verify_g13.py checks it), so the largest class of a 5-colouring has 34,
    35 or 36 points;
  - the lemma of notes/g13_chi.md, section 2 (the docstring of scripts/g13/chi/plan_C.py): if G_13 has a 5-colouring
    whose largest class, maximised over all 5-colourings, has s points, then formula F36, F35 or F34 has a solution,
    for s = 36, 35, 34;
  - the three formulas are unsatisfiable: every leaf of their cube trees has a drat-trim VERIFIED line in the logs
    with the SHA-256 of its formula, or was split again into a tree whose leaves all are, recursively. The run checked
    the five formulas of alpha(G_13) <= 36 in the same way (the case s >= 37).

This script checks the first fact and the last. It checks the archive (by default
certificates/g13_chi_certlogs.tar.gz) against its SHA256SUMS file (g13_chi_SHA256SUMS.txt beside it): the SHA-256 of
the archive, and of every file in it. It unpacks the archive into a temporary copy of scripts/g13/chi, and runs there
  - `share_driver.py regen`, which writes the eight case formulas with the code, checks them against cases/SHA256SUMS,
    and writes every re-split formula from its parent formula and cube;
  - `verify_plan_D.py`, which writes the case formulas with the code again and compares them, checks every cube file
    for a cover, and accepts a leaf only with a VERIFIED line for the SHA-256 of its formula, or a re-split that is
    fully certified.
The exit status is 0 only if verify_plan_D.py ends with PLAN D FULLY CERTIFIED. No solver runs: the logs record that
drat-trim verified a DRAT proof of each leaf, and the proofs were not kept (notes/g13_chi.md, section 7, says how to
refute a leaf again). A log line in which kissat found a leaf satisfiable is a failure too. Each re-split formula takes about 9 MB on the disk while the check runs; the script stops at once
if the disk has too little room (--force goes on anyway).

With --stats, it also prints the numbers that notes/g13_chi.md reports: for each case, the leaves of its cube tree
and of the re-splits, level by level, and the certified leaves; the log lines, the timeouts, and the kissat and
drat-trim time that the logs record. It follows the rule of verify_plan_D.py (a VERIFIED line first, then a re-split
CASE_leafI, then CASE_splitI), computing the SHA-256 of every leaf formula in memory, and reports a difference from
the counts that verify_plan_D.py prints. With --fill NOTE, it writes the numbers into the placeholders {{NAME}} of
NOTE and lists the placeholders left. With --no-verify, it prints the numbers without running regen and
verify_plan_D.py, so no re-split formula is written; the exit status then says only whether the tree it counted is
complete.

usage: python3 scripts/verify_g13_chi.py [ARCHIVE] [--sums FILE] [--workdir DIR] [--keep] [--force]
                                         [--stats] [--fill NOTE] [--no-verify]
       python3 scripts/verify_g13_chi.py --dir ASSEMBLY [--stats] [--fill NOTE] [--no-verify]
With --dir, ASSEMBLY is a merged copy of scripts/g13/chi (README.md there, "Checking the result"); regen and
verify_plan_D.py run in it.
"""
import argparse
import collections
import glob
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CHI = os.path.join(HERE, "g13", "chi")
ARCHIVE = os.path.join(ROOT, "certificates", "g13_chi_certlogs.tar.gz")
SUMS_NAME = "g13_chi_SHA256SUMS.txt"
COLOURINGS = os.path.join(ROOT, "data", "small_plane_colourings.json")
TREES = ("E37_B", "F36", "F35", "F34")                 # the case formulas with a cube tree, in verify_plan_D.py's order
ROSETTES = (6, 7, 9, 11)                               # the formulas E37_A_c, refuted without a cube tree
CASE = r"(?:E37_B|F36|F35|F34)"
MEMBERS = [re.compile(rf"cases/{CASE}(?:_(?:leaf|split)\d+)*\.certlog"),
           re.compile(rf"cases/{CASE}(?:_(?:leaf|split)\d+)+\.icnf"),
           re.compile(r"logs/E37_A\.log")]
LINE = re.compile(r"(\S+): sha256 ([0-9a-f]{64}); kissat (UNSAT|SAT|UNKNOWN) in ([\d.]+) s \(proof ([\d.]+) MB\); "
                  r"drat-trim (VERIFIED|NOT VERIFIED|not run)(?:, \d+ of \d+ lemmas in core)? in ([\d.]+) s$")
LIMITS = (120, 1200, 3600)                             # the kissat limits of share_driver.py
FORMULA_BYTES = 9.2e6                                  # the largest case formula, F36, has 9 199 090 bytes


def say(line):
    print(line, flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_colouring():
    """the stored 6-colouring of G_13, on edges rebuilt here: (True or False, report)"""
    col = json.load(open(COLOURINGS))["G_q"]["13"]
    edges = [(v, w) for v in range(169) for w in range(v + 1, 169)
             if ((v // 13 - w // 13) ** 2 - 2 * (v % 13 - w % 13) ** 2) % 13 == 1]
    bad = sum(col[v] == col[w] for v, w in edges)
    ok = len(col) == 169 and set(col) == set(range(6)) and len(edges) == 1183 and bad == 0
    return ok, (f"the stored colouring of G_13 (data/small_plane_colourings.json): {len(set(col))} colours, "
                f"{bad} of the {len(edges)} edges monochromatic")


# ---------------------------------------------------------------------------------------------------------- archive

def read_sums(path):
    """name -> SHA-256, from a file in the format of sha256sum"""
    out = {}
    for line in open(path):
        m = re.fullmatch(r"([0-9a-f]{64})  (\S+)\n?", line)
        if not m or m.group(2) in out:
            raise ValueError(f"{path}: unexpected line {line!r}")
        out[m.group(2)] = m.group(1)
    return out


def unpack(archive, sums, work, failures):
    """write the files of the archive into work, after checking every name, and the SHA-256 of the archive and of
    each file against sums (if given); returns the number of files"""
    if sums is not None:
        want = sums.get(os.path.basename(archive))
        got = sha256_file(archive)
        if want != got:
            failures.append(f"{archive}: sha256 {got}, not the {want} of its SHA256SUMS file")
            return 0
        say(f"{archive}: sha256 {got}, as in its SHA256SUMS file")
    n = 0
    seen = set()
    with tarfile.open(archive, "r:gz") as tar:
        for info in tar:
            name = info.name
            if not info.isfile() or not any(p.fullmatch(name) for p in MEMBERS) or name in seen:
                failures.append(f"{archive}: unexpected member {name!r}")
                return n
            seen.add(name)
            data = tar.extractfile(info).read()
            if sums is not None and sums.get(name) != hashlib.sha256(data).hexdigest():
                failures.append(f"{archive}: {name} does not have the SHA-256 of the SHA256SUMS file")
                return n
            with open(os.path.join(work, name), "wb") as fh:
                fh.write(data)
            n += 1
    if sums is not None:
        missing = sorted(set(sums) - seen - {os.path.basename(archive)})
        if missing:
            failures.append(f"{archive}: {len(missing)} files of the SHA256SUMS file are missing, e.g. {missing[0]}")
    return n


def copy_code(chi, work):
    """a copy of scripts/g13/chi without anything a run writes (formulas, logs, re-splits, proofs, flags)"""
    def ignore(folder, names):
        rel = os.path.relpath(folder, chi)
        out = {n for n in names if n in ("__pycache__", "proofs", "logs", "progress.txt", "DONE", "SAT", "FAILED")
               or n.endswith((".cnf", ".cnf.tmp", ".certlog", ".todo", ".model", ".tar.gz"))}
        if rel == "cases":
            out |= {n for n in names if re.fullmatch(rf"{CASE}(?:_(?:leaf|split)\d+)+\.icnf", n)}
        return out
    shutil.copytree(chi, work, ignore=ignore)
    os.makedirs(os.path.join(work, "logs"), exist_ok=True)


def room_for_regen(work):
    """(bytes that share_driver.py regen will write, bytes free)"""
    need = 0
    for icnf in glob.glob(os.path.join(work, "cases", "*_leaf*.icnf")) + glob.glob(os.path.join(work, "cases", "*_split*.icnf")):
        if not os.path.exists(icnf[:-5] + ".cnf"):
            need += FORMULA_BYTES
    for name in ("E37_A6", "E37_A7", "E37_A9", "E37_A11") + TREES:
        if not os.path.exists(os.path.join(work, "cases", name + ".cnf")):
            need += FORMULA_BYTES
    return int(need), shutil.disk_usage(work).free


# ------------------------------------------------------------------------------------------------------- statistics

def parse_log(path):
    """the lines of a log of certify.py, parsed: (name, sha256, kissat verdict, kissat s, proof MB, drat-trim verdict,
    drat-trim s); and the lines of another form"""
    lines, odd = [], []
    for text in (open(path) if os.path.exists(path) else []):
        m = LINE.match(text.rstrip("\n"))
        if m:
            g = m.groups()
            lines.append((g[0], g[1], g[2], float(g[3]), float(g[4]), g[5], float(g[6])))
        elif text.strip():
            odd.append(text)
    return lines, odd


def cubes_of(path):
    """the leaves of a cube file, as the strings of their literals (as certify.py writes them), and how many are
    'c closed'"""
    out, closed = [], 0
    for line in open(path):
        if line.startswith("a "):
            out.append(line.split()[1:-1])
        elif line.startswith("c closed"):
            out.append(line.split()[2:-1])
            closed += 1
    return out, closed


class Tree:
    """what a certified (sub)tree contains: per level, the formulas, leaves, closed leaves, leaves with a VERIFIED
    line, and leaves split again (by kind); the VERIFIED lines it uses; its re-split formulas"""

    def __init__(self):
        self.levels = collections.defaultdict(collections.Counter)
        self.used, self.resplits = [], []

    def add(self, other):
        for k, c in other.levels.items():
            self.levels[k].update(c)
        self.used += other.used
        self.resplits += other.resplits


class Walker:
    """the certified tree of each case formula, by the rule of verify_plan_D.py, with the SHA-256 of every leaf
    formula computed in memory: the formula of a re-split is its parent's leaf formula, so a leaf of it is the case
    formula with the cubes of the path to it as unit clauses"""

    def __init__(self, work, check_cover):
        self.work, self.check_cover, self.logs, self.why = work, check_cover, {}, []

    def verified(self, stem):
        """leaf name -> SHA-256 -> the first VERIFIED line of the log of stem"""
        if stem not in self.logs:
            ok = collections.defaultdict(dict)
            for line in parse_log(os.path.join(self.work, stem + ".certlog"))[0]:
                if line[2] == "UNSAT" and line[5] == "VERIFIED":
                    ok[line[0]].setdefault(line[1], line)
            self.logs[stem] = ok
        return self.logs[stem]

    def walk(self, stem, nv, nc, body, extra, level):
        """the Tree under the cube file stem.icnf of the formula 'p cnf nv nc' + body + extra, or None"""
        icnf = os.path.join(self.work, stem + ".icnf")
        if not os.path.exists(icnf):
            self.why.append(f"{stem}.icnf is missing")
            return None
        cubes, closed = cubes_of(icnf)
        if not self.check_cover([tuple(map(int, c)) for c in cubes]):
            self.why.append(f"{stem}.icnf is not a cover")
            return None
        tree, ok, tag, prefix = Tree(), self.verified(stem), os.path.basename(stem), {}
        here = tree.levels[level]
        here.update(formulas=1, leaves=len(cubes), closed=closed)
        for i, cube in enumerate(cubes):
            units = "".join(f"{lit} 0\n" for lit in cube).encode()
            if len(cube) not in prefix:
                h = hashlib.sha256(f"p cnf {nv} {nc + len(cube)}\n".encode())
                h.update(body)
                h.update(extra)
                prefix[len(cube)] = h
            h = prefix[len(cube)].copy()
            h.update(units)
            line = ok.get(f"{tag}_leaf{i}", {}).get(h.hexdigest())
            if line:
                here["verified"] += 1
                tree.used.append(line)
                continue
            sub = None
            for kind in ("leaf", "split"):                      # verify_plan_D.py tries them in this order
                name = f"{stem}_{kind}{i}"
                if sub is None and os.path.exists(os.path.join(self.work, name + ".icnf")):
                    sub = self.walk(name, nv, nc + len(cube), body, extra + units, level + 1)
                    if sub is not None:
                        here[f"re-split ({kind})"] += 1
                        tree.resplits.append(name)
            if sub is None:
                self.why.append(f"{tag}_leaf{i}: not VERIFIED and not re-split into a certified case")
                return None
            here["re-split"] += 1
            tree.add(sub)
        return tree


def statistics(work):
    """the numbers of notes/g13_chi.md from the merged logs in work (the case formulas must be in work/cases):
    (a report, a dict of placeholders, True if every case is fully certified)"""
    spec = importlib.util.spec_from_file_location("g13chi_cuber2_stats", os.path.join(work, "cuber2.py"))
    cuber2 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cuber2)
    sys.setrecursionlimit(100000)
    walker = Walker(work, cuber2.check_cover)
    sums = dict(reversed(line.split()) for line in open(os.path.join(work, "cases", "SHA256SUMS")) if line.strip())
    report, ph, complete = [], {}, True

    a_lines = parse_log(os.path.join(work, "logs", "E37_A.log"))[0]
    a_ok = {line[1] for line in a_lines if line[2] == "UNSAT" and line[5] == "VERIFIED"}
    a_used = [next(line for line in a_lines if line[1] == sums[f"E37_A{c}.cnf"] and line[5] == "VERIFIED")
              for c in ROSETTES if sums[f"E37_A{c}.cnf"] in a_ok]
    complete &= len(a_used) == 4
    report.append(f"E37_A6, E37_A7, E37_A9, E37_A11: {len(a_used)} of 4 with a VERIFIED line (logs/E37_A.log)")

    total, used, trees = Tree(), list(a_used), {}
    for case in TREES:
        path = os.path.join(work, "cases", case + ".cnf")
        with open(path, "rb") as fh:
            head, body = fh.read().split(b"\n", 1)
        if hashlib.sha256(head + b"\n" + body).hexdigest() != sums[case + ".cnf"]:
            raise SystemExit(f"{path} does not have the SHA-256 of cases/SHA256SUMS")
        nv, nc = map(int, head.split()[2:4])
        tree = walker.walk(f"cases/{case}", nv, nc, body, b"", 0)
        if tree is None:
            complete = False
            report.append(f"{case}: NOT fully certified: {walker.why[-1]}")
            continue
        trees[case] = tree
        total.add(tree)
        top = tree.levels[0]
        certified = sum(c["verified"] for c in tree.levels.values())
        depth = max(tree.levels)
        report.append(f"{case}: {top['leaves']} leaves ({top['closed']} closed by unit propagation); "
                      f"{top['verified']} VERIFIED, {top['re-split']} split again ({top['re-split (split)']} as "
                      f"CASE_splitI, {top['re-split (leaf)']} as CASE_leafI); certified in all: {certified}"
                      + (f"; the deepest re-split is at level {depth}" if depth else ""))
        for level in sorted(tree.levels)[1:]:
            c = tree.levels[level]
            report.append(f"    level {level}: {c['formulas']} re-split formulas, {c['leaves']} leaves "
                          f"({c['closed']} closed); {c['verified']} VERIFIED, {c['re-split']} split again")
        ph.update({f"{case}_TOP_VERIFIED": top["verified"], f"{case}_TOP_RESPLIT": top["re-split"],
                   f"{case}_CERTIFIED": certified, f"{case}_DEPTH": depth,
                   f"{case}_RESPLITS": sum(c["formulas"] for k, c in tree.levels.items() if k)})

    used += total.used
    alllines, odd = [], []
    for log in sorted(glob.glob(os.path.join(work, "cases", "*.certlog"))) + [os.path.join(work, "logs", "E37_A.log")]:
        lines, other = parse_log(log)
        alllines += lines
        odd += other
    verified = sum(1 for x in alllines if x[2] == "UNSAT" and x[5] == "VERIFIED")
    unknown = [x for x in alllines if x[2] == "UNKNOWN"]
    by_limit = collections.Counter(next((t for t in LIMITS if x[3] >= 0.95 * t and x[3] < 1.5 * t), "other")
                                   for x in unknown)
    sat = sum(1 for x in alllines if x[2] == "SAT")
    ph["SAT_LINES"] = sat
    not_verified = sum(1 for x in alllines if x[2] == "UNSAT" and x[5] != "VERIFIED")
    files = set(glob.glob(os.path.join(work, "cases", "*_leaf*.icnf")) + glob.glob(os.path.join(work, "cases", "*_split*.icnf")))
    unused = len(files) - len(total.resplits)
    hours = lambda xs, k: sum(x[k] for x in xs) / 3600
    certified = len(used)
    if complete:
        report.append(f"all cases: {certified} leaves and formulas certified, each by one VERIFIED line (4 formulas "
                      f"E37_A_c and {certified - 4} leaves); {len(total.resplits)} re-split formulas")
        report.append(f"the lines that certify them: kissat {group(round(3600 * hours(used, 3)))} s "
                      f"({hours(used, 3):.1f} h), drat-trim {group(round(3600 * hours(used, 6)))} s "
                      f"({hours(used, 6):.1f} h); "
                      f"the longest kissat run {max(x[3] for x in used):.1f} s; DRAT proofs of "
                      f"{sum(x[4] for x in used) / 1e3:.1f} GB in all, deleted after the check")
    report.append(f"all logs: {len(alllines)} lines of certify.py ({len(odd)} other lines); {verified} VERIFIED, "
                  f"{len(unknown)} timeouts ({', '.join(f'{by_limit[t]} at {t} s' for t in LIMITS)}"
                  + (f", {by_limit['other']} stopped early" if by_limit["other"] else "") + f"), {sat} SAT, "
                  f"{not_verified} UNSAT but not VERIFIED; kissat {group(round(3600 * hours(alllines, 3)))} s "
                  f"({hours(alllines, 3):.1f} h), drat-trim {group(round(3600 * hours(alllines, 6)))} s "
                  f"({hours(alllines, 6):.1f} h)")
    report.append(f"cube files of re-splits: {len(files)}, of which {unused} are not used by the certified trees")
    plan_d = [cubes_of(f)[1:] + (len(cubes_of(f)[0]),) for f in files if re.search(r"_leaf\d+\.icnf$", f)]
    report.append(f"re-splits of plan D (CASE_leafI): {len(plan_d)}, of which {plan_d.count((8, 9))} have 9 leaves, "
                  f"8 of them closed")
    if not complete:
        return report, ph, False

    rows = []
    for case in ("F36", "F35", "F34"):
        for level in sorted(trees[case].levels)[1:]:
            c = trees[case].levels[level]
            rows.append(f"| `{case}` | {level} | {group(c['formulas'])} | {group(c['leaves'])} | "
                        f"{group(c['verified'])} | {group(c['re-split'])} |")
    colouring = [ph[f"F3{s}_CERTIFIED"] for s in (6, 5, 4)]
    ph.update({"TOTAL_CERTIFIED": certified, "DRAT_PROOFS": certified, "COLOURING_CERTIFIED": sum(colouring),
               "RESPLITS": len(total.resplits), "MAX_DEPTH": max(max(t.levels) for t in trees.values()),
               "KISSAT_HOURS": f"{hours(used, 3):,.1f}".replace(",", " "),
               "DRAT_TRIM_HOURS": f"{hours(used, 6):,.1f}".replace(",", " "),
               "MAX_KISSAT_S": round(max(x[3] for x in used)), "PROOF_GB": round(sum(x[4] for x in used) / 1e3),
               "ALL_LINES": len(alllines), "ALL_VERIFIED": verified, "TIMEOUTS": len(unknown),
               "TIMEOUTS_120": by_limit[120], "TIMEOUTS_1200": by_limit[1200], "TIMEOUTS_3600": by_limit[3600],
               "ALL_KISSAT_HOURS": round(hours(alllines, 3)), "ALL_DRAT_TRIM_HOURS": round(hours(alllines, 6)),
               "UNUSED_RESPLITS": unused, "REGEN_GB": f"{len(files) * FORMULA_BYTES / 1e9:.1f}",
               "TREE_TABLE": "\n".join(rows) if rows else "| (none) | | | | | |"})
    return report, ph, True


def group(n):
    """1234567 -> '1 234 567', as the notes write numbers"""
    return f"{n:,}".replace(",", " ")


def fill(path, ph):
    """write the placeholders {{NAME}} of ph into the file; returns the placeholders left in it"""
    text = open(path).read()
    for k, v in ph.items():
        text = text.replace("{{" + k + "}}", group(v) if isinstance(v, int) else str(v))
    with open(path, "w") as fh:
        fh.write(text)
    return sorted(set(re.findall(r"\{\{([A-Z0-9_]+)\}\}", text)))


# ------------------------------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description="Recheck chi(G_13) = 6 from the certificate archive.")
    ap.add_argument("archive", nargs="?", default=ARCHIVE)
    ap.add_argument("--sums", default=None, help=f"default: {SUMS_NAME} beside the archive")
    ap.add_argument("--dir", default=None, help="a merged assembly to check in place, instead of an archive")
    ap.add_argument("--chi", default=CHI, help="the code to copy (default scripts/g13/chi)")
    ap.add_argument("--workdir", default=None, help="where to make the temporary copy (default: the system's)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary copy")
    ap.add_argument("--force", action="store_true", help="go on even if the disk seems to have too little room")
    ap.add_argument("--stats", action="store_true", help="print the numbers of notes/g13_chi.md")
    ap.add_argument("--fill", default=None, metavar="NOTE", help="write the numbers into the placeholders of NOTE")
    ap.add_argument("--no-verify", action="store_true", help="only the numbers: no regen, no verify_plan_D.py")
    a = ap.parse_args()
    a.stats = a.stats or bool(a.fill)
    if a.no_verify and not a.stats:
        ap.error("--no-verify needs --stats or --fill")
    failures = []
    ok, report = check_colouring()
    say(report)
    if not ok:
        failures.append(report)

    tmp = None
    if a.dir:
        work = os.path.abspath(a.dir)
    else:
        sums_path = a.sums or os.path.join(os.path.dirname(os.path.abspath(a.archive)), SUMS_NAME)
        if not os.path.exists(a.archive):
            say(f"NOT CONFIRMED: there is no archive {a.archive} (FINALIZE.md says how it is made)")
            sys.exit(1)
        sums = read_sums(sums_path) if os.path.exists(sums_path) else None
        if sums is None:
            say(f"no {sums_path}: the archive is not checked against a SHA256SUMS file")
        tmp = tempfile.mkdtemp(prefix="g13_chi_", dir=a.workdir)
        work = os.path.join(tmp, "chi")
        copy_code(os.path.abspath(a.chi), work)
        n = unpack(a.archive, sums, work, failures)
        say(f"{n} files unpacked into a copy of {a.chi}: {work}")
    try:
        counts, complete = None, None
        if not failures and not a.no_verify:
            counts = verify(work, a.force, failures)
        if a.stats and (not failures or counts is not None):     # after a failed verify_plan_D.py too, to see why
            complete = report_statistics(a, work, counts or {}, failures)
    finally:
        if tmp and not a.keep:
            shutil.rmtree(tmp, ignore_errors=True)
        elif tmp:
            say(f"kept {work}")

    if failures:
        say(f"NOT CONFIRMED: {len(failures)} failures, the first: {failures[0]}")
        sys.exit(1)
    if a.no_verify:
        say("STATISTICS ONLY: regen and verify_plan_D.py were not run; the trees counted are "
            + ("complete." if complete else "NOT complete."))
        sys.exit(0 if complete else 1)
    say("CONFIRMED: the stored 6-colouring, the eight case formulas against the code, every cover, and a drat-trim "
        "VERIFIED line for every leaf of every cube tree (PLAN D FULLY CERTIFIED). With alpha(G_13) = 36 "
        "(scripts/verify_g13.py) and the lemma of notes/g13_chi.md, chi(G_13) = 6.")


def verify(work, force, failures):
    """share_driver.py regen and verify_plan_D.py in work: the counts of certified leaves that verify_plan_D.py prints
    (None if it did not run); a failure goes to failures"""
    need, free = room_for_regen(work)
    say(f"share_driver.py regen will write about {need / 1e9:.1f} GB; {free / 1e9:.1f} GB are free")
    if need + 1e8 > free and not force:
        failures.append("too little room on the disk for the re-split formulas (give --workdir, or --force)")
        return None
    r = subprocess.run([sys.executable, "share_driver.py", "regen"], cwd=work, capture_output=True, text=True)
    say("share_driver.py regen: " + (r.stdout + r.stderr).strip().replace("\n", "\n    "))
    if r.returncode:
        failures.append(f"share_driver.py regen: exit status {r.returncode}")
        return None
    r = subprocess.run([sys.executable, "verify_plan_D.py"], cwd=work, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    say("verify_plan_D.py:\n    " + out.replace("\n", "\n    "))
    if r.returncode or out.splitlines()[-1:] != ["PLAN D FULLY CERTIFIED"]:
        failures.append(f"verify_plan_D.py: exit status {r.returncode}, last line {out.splitlines()[-1:]}")
    sat = [line for log in glob.glob(os.path.join(work, "cases", "*.certlog")) + [os.path.join(work, "logs", "E37_A.log")]
           for line in parse_log(log)[0] if line[2] == "SAT"]
    if sat:
        failures.append(f"kissat found {len(sat)} formulas satisfiable, e.g. {sat[0][0]}: check the model "
                        f"(check_colouring.py)")
    return {m.group(1): int(m.group(2)) for m in re.finditer(r"^(\S+) (\d+) leaves VERIFIED", r.stdout, re.M)}


def report_statistics(a, work, counts, failures):
    """print the statistics of the trees in work, compare them with the counts of verify_plan_D.py, and fill the
    note; returns whether every case is fully certified"""
    r = subprocess.run([sys.executable, "-c", "import share_driver; share_driver.formulas()"], cwd=work,
                       capture_output=True, text=True)          # the case formulas, written if missing, and checked
    if r.returncode:
        failures.append(f"the case formulas: {(r.stdout + r.stderr).strip()}")
        return False
    lines, ph, complete = statistics(work)
    if ph.get("SAT_LINES") and a.no_verify:
        failures.append(f"kissat found {ph['SAT_LINES']} formulas satisfiable: check the models (check_colouring.py)")
    if a.dir is None and os.path.exists(a.archive):
        ph.update({"ARCHIVE_SHA256": sha256_file(a.archive), "ARCHIVE_MB": f"{os.path.getsize(a.archive) / 1e6:.1f}",
                   "ARCHIVE_FILES": sum(1 for _ in tarfile.open(a.archive, "r:gz"))})
    say("statistics (notes/g13_chi.md, section 5):\n    " + "\n    ".join(lines))
    say("placeholders:\n" + "\n".join(f"    {{{{{k}}}}} {group(v) if isinstance(v, int) else v}".replace("\n", "\n        ")
                                       for k, v in sorted(ph.items())))
    for case, n in counts.items():
        if ph.get(f"{case}_CERTIFIED") != n:
            failures.append(f"{case}: verify_plan_D.py counts {n} leaves, the statistics {ph.get(f'{case}_CERTIFIED')}")
    if a.fill:
        left = fill(a.fill, ph)
        say(f"{a.fill}: placeholders written; left: {', '.join(left) if left else 'none'}")
    return complete


if __name__ == "__main__":
    main()
