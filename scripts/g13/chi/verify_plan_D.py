"""Check that plan D (run_plan_D.sh) is fully certified: exit 0 only if every leaf of every cover is certified.

For each case formula it recomputes the cover of its cube file (cuber2.check_cover), and for every leaf the SHA-256
of the leaf formula that certify.py writes; the leaf counts only if a drat-trim VERIFIED line in the case's log has
its name and that SHA-256, or if it was re-split (resplit.sh: CASE_leafI.cnf; plan E: CASE_splitI.cnf; either
byte-identical to the leaf formula) and one such re-split case is itself fully certified, recursively. It also regenerates every case formula with the code
and compares it with the one that was certified, and checks the four Part A formulas of E37 against logs/E37_A.log.
Written after an independent review of plan D (review_plan_D.md, section 4), which it follows.
usage (in g13/, after run_plan_D.sh): python3 verify_plan_D.py
"""
import hashlib
import os
import re
import subprocess
import sys
import tempfile

sys.setrecursionlimit(100000)
from cuber2 import check_cover  # noqa: E402


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def verified(log):
    """name -> set of SHA-256 with a drat-trim VERIFIED line"""
    ok = {}
    for ln in (open(log) if os.path.exists(log) else []):
        m = re.match(r"(\S+): sha256 ([0-9a-f]{64}); .*; drat-trim VERIFIED", ln)
        if m:
            ok.setdefault(m.group(1), set()).add(m.group(2))
    return ok


def cubes_of(icnf):
    out = []
    for line in open(icnf):
        if line.startswith("a "):
            out.append(line.split()[1:-1])
        elif line.startswith("c closed"):
            out.append(line.split()[2:-1])
    return out


FAILS = []


def check_case(cnf):
    """the number of VERIFIED leaves under the cube file of cnf, re-splits included; None (and a reason in FAILS)
    when some leaf is neither VERIFIED nor re-split into a fully certified case"""
    stem = cnf[:-4]
    tag = os.path.basename(stem)
    head, body = open(cnf, "rb").read().split(b"\n", 1)
    nv, nc = map(int, head.split()[2:4])
    if not os.path.exists(stem + ".icnf"):
        FAILS.append(f"{stem}.icnf is missing")
        return None
    cubes = cubes_of(stem + ".icnf")
    if not check_cover([tuple(map(int, c)) for c in cubes]):
        FAILS.append(f"{stem}.icnf is not a cover")
        return None
    ok, pre, n = verified(stem + ".certlog"), {}, 0
    for i, cube in enumerate(cubes):
        if len(cube) not in pre:
            pre[len(cube)] = hashlib.sha256(f"p cnf {nv} {nc + len(cube)}\n".encode() + body)
        h = pre[len(cube)].copy()
        h.update("".join(f"{lit} 0\n" for lit in cube).encode())
        name, H = f"{tag}_leaf{i}", h.hexdigest()          # exactly the leaf formula certify.py hashes
        if H in ok.get(name, ()):
            n += 1
            continue
        m = None
        for sub in (f"{stem}_leaf{i}.cnf", f"{stem}_split{i}.cnf"):   # resplit.sh, plan E: the leaf formula itself
            if m is None and os.path.exists(sub) and sha(sub) == H:
                m = check_case(sub)
        if m is None:
            FAILS.append(f"{name}: not VERIFIED and not re-split into a certified case")
            return None
        n += m
    return n


GEN = {f"E37_A{c}": f"enum_cert.py {{}} 37 --L 0 --rosette {c}" for c in (6, 7, 9, 11)}
GEN.update({"E37_B": "enum_cert.py {} 37 --norosette",
            "F36": "g13cnf.py {} --big0 36 --cap 36 --dom0 --lex0 25 --vp 1234 --domcap 36",
            "F35": "g13cnf.py {} --big0 35 --cap 35 --dom0 --lex0 25 --vp 1234 --domcap 35",
            "F34": "g13cnf.py {} --rigid34 --lex0 25"})
with tempfile.TemporaryDirectory() as d:                  # the certified formulas are what the code writes
    for k, cmd in GEN.items():
        out = os.path.join(d, k + ".cnf")
        subprocess.run(["python3"] + cmd.format(out).split(), check=True, capture_output=True)
        if sha(out) != sha(f"cases/{k}.cnf"):
            sys.exit(f"FAIL cases/{k}.cnf differs from the code's output")
okA = {h for s in verified("logs/E37_A.log").values() for h in s}
for c in (6, 7, 9, 11):
    if sha(f"cases/E37_A{c}.cnf") not in okA:
        sys.exit(f"FAIL E37_A{c} not VERIFIED")
for case in ("E37_B", "F36", "F35", "F34"):
    n = check_case(f"cases/{case}.cnf")
    if n is None:
        sys.exit(f"FAIL {case}: {FAILS[-1]}")
    print(case, n, "leaves VERIFIED (re-split sub-leaves included)")
print("PLAN D FULLY CERTIFIED")
