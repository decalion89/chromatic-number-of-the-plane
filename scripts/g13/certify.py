"""Certify formulas as unsatisfiable: kissat writes a binary DRAT proof, drat-trim checks it against the formula,
then the proof (and any temporary leaf formula) is deleted.  One log line per formula: name, sha256 of the exact
formula text, kissat verdict and time, proof size, drat-trim verdict (VERIFIED only if drat-trim prints
's VERIFIED'), lemmas in core, drat-trim time.  A SAT answer is logged and kissat's model saved as NAME.model.

usage: certify.py --log LOG [--jobs N] [--time T] FILE.cnf ...
       certify.py --log LOG [--jobs N] [--time T] --base CASE.cnf --cubes CUBES.icnf [--only i,j,...]
In cube mode, leaf i (line i of the cube file, 'a ...' or 'c closed ...') is CASE.cnf plus its literals as unit
clauses; the leaf file is written to --workdir, certified and deleted.  Check the cover first (cuber2.py does,
and check_cover() in cuber2.py re-checks a file).  Every solver call runs under 'nice -n 19'.
"""
import argparse
import hashlib
import os
import re
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
KISSAT = os.environ.get("KISSAT", "kissat")              # kissat 4.0.4 and drat-trim 2e3b2dc (worker_setup.sh)
DRAT_TRIM = os.environ.get("DRAT_TRIM", "drat-trim")


def certify(path, name, time_limit, workdir, drat_timeout=200000):
    text = open(path, "rb").read()
    sha = hashlib.sha256(text).hexdigest()
    stem = os.path.join(workdir, name)
    proof = stem + ".drat"
    t0 = time.time()
    r = subprocess.run(["nice", "-n", "19", KISSAT, "--unsat", f"--time={time_limit}", path, proof],
                       capture_output=True, text=True)
    t1 = time.time()
    verdict = {10: "SAT", 20: "UNSAT"}.get(r.returncode, "UNKNOWN")
    checked, core, psize = "not run", "", 0
    if verdict == "SAT":
        with open(stem + ".model", "w") as fh:
            fh.write(r.stdout)
    if os.path.exists(proof):
        psize = os.path.getsize(proof)
    if verdict == "UNSAT":
        d = subprocess.run(["nice", "-n", "19", DRAT_TRIM, path, proof, "-t", str(drat_timeout)],
                           capture_output=True, text=True)
        out = d.stdout + d.stderr
        checked = "VERIFIED" if re.search(r"^s VERIFIED", out, re.M) else "NOT VERIFIED"
        m = re.search(r"(\d+) of (\d+) lemmas in core", out)
        core = f", {m.group(1)} of {m.group(2)} lemmas in core" if m else ""
    t2 = time.time()
    if os.path.exists(proof):
        os.remove(proof)
    return (f"{name}: sha256 {sha}; kissat {verdict} in {t1 - t0:.1f} s (proof {psize / 1e6:.1f} MB); "
            f"drat-trim {checked}{core} in {t2 - t1:.1f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*")
    ap.add_argument("--log", required=True)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--time", type=int, default=1200, help="kissat time limit per formula (s)")
    ap.add_argument("--workdir", default=os.path.join(HERE, "proofs"))
    ap.add_argument("--base", default=None)
    ap.add_argument("--cubes", default=None)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)
    lock = threading.Lock()

    def log(line):
        with lock:
            print(line, flush=True)
            with open(a.log, "a") as fh:
                fh.write(line + "\n")

    if a.cubes:
        head, body = open(a.base).read().split("\n", 1)
        nv, nc = map(int, head.split()[2:4])
        cubes = []
        for line in open(a.cubes):
            if line.startswith("a "):
                cubes.append(line.split()[1:-1])
            elif line.startswith("c closed"):
                cubes.append(line.split()[2:-1])
        todo = [int(x) for x in a.only.split(",") if x] or range(len(cubes))
        tag = os.path.splitext(os.path.basename(a.cubes))[0]
        done = set()
        if os.path.exists(a.log):          # resumable: skip a leaf only if a VERIFIED line has its name AND the
            for l in open(a.log):          # sha256 of the leaf formula about to be written (not a stale log)
                m = re.match(r"(\S+): sha256 ([0-9a-f]{64}); .*drat-trim VERIFIED", l)
                if m:
                    done.add(m.groups())

        def run(i):
            name = f"{tag}_leaf{i}"
            text = f"p cnf {nv} {nc + len(cubes[i])}\n" + body + "".join(f"{l} 0\n" for l in cubes[i])
            if (name, hashlib.sha256(text.encode()).hexdigest()) in done:
                return
            path = os.path.join(a.workdir, name + ".cnf")
            with open(path, "w") as fh:
                fh.write(text)
            try:
                log(certify(path, name, a.time, a.workdir))
            finally:
                os.remove(path)
        with ThreadPoolExecutor(a.jobs) as ex:
            list(ex.map(run, todo))
    else:
        with ThreadPoolExecutor(a.jobs) as ex:
            list(ex.map(lambda p: log(certify(p, os.path.basename(p), a.time, a.workdir)), a.files))


if __name__ == "__main__":
    main()
