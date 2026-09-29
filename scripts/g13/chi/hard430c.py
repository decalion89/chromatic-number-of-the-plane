"""Certify level-1 sub-leaves of F34 leaf 430 of G13 on any machine, exactly as hard430b.py does on the coordinating
session's machine: re-split style 'split', cuts on colour 0 first, then colours 1 and 2 when colour 0 is exhausted
(work/vars_c012_lex.txt), kissat limits [120] + [1200] * 10. cases/F34_split430.icnf and .certlog are the level-1
split of leaf 430 made by share 33 (claude/g13-share-33).

usage (in g13-plan-d/):  python3 hard430c.py JOBS I1,I2,...
KISSAT and DRAT_TRIM (environment) as for share_driver.py. At the end it writes DONE, or SAT / FAILED as
share_driver.py does. Resumable: rerun it with the same arguments."""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)
import share_driver as sd  # noqa: E402

SHA430 = "38b9ce26685c5f9b04bf6bee95cade45d25dde94963c46b06e553df3edfb122c"   # sha256 of cases/F34_split430.cnf

jobs, idxs = int(sys.argv[1]), [int(x) for x in sys.argv[2].split(",")]
sd.formulas()
for f in ("DONE", "SAT", "FAILED"):
    if os.path.exists(f):
        os.remove(f)
os.makedirs("logs", exist_ok=True)
cnf = "cases/F34_split430.cnf"
text = sd.leaf_text("cases/F34.cnf", sd.cubes_of("cases/F34.icnf")[430])
assert hashlib.sha256(text.encode()).hexdigest() == SHA430, "the leaf formula of F34 leaf 430 differs"
if not (os.path.exists(cnf) and sd.sha_file(cnf) == SHA430):
    with open(cnf + ".tmp", "w") as fh:
        fh.write(text)
    os.replace(cnf + ".tmp", cnf)
sd.VARS["F34"] = "work/vars_c012_lex.txt"
sd.say(f"F34_split430 leaves {idxs}: start, {jobs} jobs, cuts on colour 0, then 1 and 2")
sd.resolve("F34", cnf, idxs, 1, jobs, style="split", limits=[120] + [1200] * 10)
sd.finish("DONE", f"F34_split430 leaves {idxs} certified (re-splits included)")
