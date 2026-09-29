"""Help one share of plan E with the level-1 sub-leaves of one F34 leaf, from the other end: for F34 leaf LEAF of
share S, certify the given sub-leaves I of F34_splitLEAF in the order given, exactly as the plan-E branch of
share_driver.py does (re-split style 'split', LIMITS_E, the share's own cut variables), and stop as soon as the share's
own branch holds the re-split of the next sub-leaf (cases/F34_splitLEAF_splitI.icnf). cases/F34_splitLEAF.icnf and
.certlog must first be copied from the share's branch claude/g13-share-S. Both runs write the same deterministic
re-split files, so verify_plan_D.py accepts either run's certificate.

usage (in g13-plan-d/):  python3 sub_help.py S LEAF JOBS I1,I2,...
KISSAT and DRAT_TRIM (environment) as for share_driver.py. Writes DONE, or SAT / FAILED as share_driver.py does.
Resumable: rerun it with the same arguments."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)
import share_driver as sd  # noqa: E402


def owner_has(share, name):
    """True if the share's own branch on origin holds cases/<name>"""
    ref = f"refs/remotes/origin/claude/g13-share-{share}"
    r = subprocess.run(["git", "fetch", "-q", "origin", f"claude/g13-share-{share}:{ref}"], capture_output=True, text=True)
    if r.returncode != 0:
        sd.say(f"sub-helper: could not fetch claude/g13-share-{share} ({r.stderr.strip()[:120]}); going on")
        return False
    return subprocess.run(["git", "cat-file", "-e", f"{ref}:g13-plan-d/cases/{name}"], capture_output=True).returncode == 0


share, leaf, jobs = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
idxs = [int(x) for x in sys.argv[4].split(",")]
sd.formulas()
for f in ("DONE", "SAT", "FAILED"):
    if os.path.exists(f):
        os.remove(f)
os.makedirs("logs", exist_ok=True)
tag = f"F34_split{leaf}"
cnf = f"cases/{tag}.cnf"
text = sd.leaf_text("cases/F34.cnf", sd.cubes_of("cases/F34.icnf")[leaf])
with open(cnf + ".tmp", "w") as fh:
    fh.write(text)
os.replace(cnf + ".tmp", cnf)
cubes = sd.cubes_of(f"cases/{tag}.icnf")
ok, sat = sd.log_state(f"cases/{tag}.certlog")
names = {n for n, _ in ok}
probe = [i for i in range(len(cubes)) if f"{tag}_leaf{i}" in names][:3]
shas = sd.leaf_shas(cnf, cubes, probe)
if not (probe and all((f"{tag}_leaf{i}", shas[i]) in ok for i in probe)):
    sd.finish("FAILED", f"{tag}: the leaf formulas differ from the share's certlog (wrong or missing cases/{tag}.*)")
sd.say(f"sub-helper of share {share}: start, {jobs} jobs, {tag} sub-leaves {idxs}")
done = []
for i in idxs:
    if owner_has(share, f"{tag}_split{i}.icnf"):
        sd.say(f"sub-helper of share {share}: the share's own run has reached {tag} sub-leaf {i}; stopping here")
        break
    sd.resolve("F34", cnf, [i], 1, jobs, style="split", limits=sd.LIMITS_E)
    done.append(i)
    sd.say(f"sub-helper of share {share}: {tag} sub-leaf {i} certified (re-splits included)")
sd.finish("DONE", f"sub-helper of share {share}: {tag} sub-leaves {done} certified")
