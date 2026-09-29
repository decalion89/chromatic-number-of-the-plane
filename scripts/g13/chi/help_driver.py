"""Help one share of plan E from the other end: certify some F34 leaves of that share in the order given (the reverse
of the share's own order), exactly as the plan-E branch of share_driver.py does (resolve() with the re-split style
'split' and the limits LIMITS_E; F35 is complete and is not run). Before each leaf the share's own branch is fetched
from origin: if it already holds the re-split of that leaf (cases/F34_splitI.icnf), the share's own run has reached
it, and the helper stops there. Both runs write the same deterministic re-split files, so verify_plan_D.py accepts
either run's certificate for a leaf.

usage (in g13-plan-d/):  python3 help_driver.py S I1,I2,...  [--jobs N]
KISSAT and DRAT_TRIM (environment) as for share_driver.py. At the end it writes DONE, or SAT / FAILED as
share_driver.py does. Resumable: rerun it with the same arguments."""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import share_driver as sd  # noqa: E402


def owner_has(share, leaf):
    """True if the share's own branch on origin holds the re-split of F34 leaf `leaf`"""
    ref = f"refs/remotes/origin/claude/g13-share-{share}"
    r = subprocess.run(["git", "fetch", "-q", "origin", f"claude/g13-share-{share}:{ref}"], cwd=HERE,
                       capture_output=True, text=True)
    if r.returncode != 0:
        sd.say(f"helper: could not fetch claude/g13-share-{share} ({r.stderr.strip()[:120]}); going on")
        return False
    r = subprocess.run(["git", "cat-file", "-e", f"{ref}:g13-plan-d/cases/F34_split{leaf}.icnf"], cwd=HERE,
                       capture_output=True)
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("share", type=int)
    ap.add_argument("leaves")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    a = ap.parse_args()
    leaves = [int(x) for x in a.leaves.split(",")]
    plan = dict((case, idxs) for case, idxs in json.load(open(sd.path("planE.json")))[str(a.share)])
    if not set(leaves) <= set(plan["F34"]):
        sys.exit(f"leaves {sorted(set(leaves) - set(plan['F34']))} are not F34 leaves of share {a.share}")
    sd.formulas()
    for f in ("DONE", "SAT", "FAILED"):
        if os.path.exists(sd.path(f)):
            os.remove(sd.path(f))
    os.makedirs(sd.path("logs"), exist_ok=True)
    sd.say(f"helper of share {a.share}: start, {a.jobs} jobs, F34 leaves {leaves}; kissat {os.environ.get('KISSAT')}, "
           f"drat-trim {os.environ.get('DRAT_TRIM')}")
    cnf = sd.path("cases", "F34.cnf")
    cubes = sd.cubes_of(sd.path("cases", "F34.icnf"))
    if not sd.check_cover([tuple(map(int, c)) for c in cubes]):
        sd.finish("FAILED", "cases/F34.icnf is not a cover")
    done = []
    for leaf in leaves:
        if owner_has(a.share, leaf):
            sd.say(f"helper of share {a.share}: the share's own run has reached F34 leaf {leaf}; stopping here")
            break
        sd.resolve("F34", cnf, [leaf], 0, a.jobs, style="split", limits=sd.LIMITS_E)
        done.append(leaf)
        sd.say(f"helper of share {a.share}: F34 leaf {leaf} certified (re-splits included)")
    sd.finish("DONE", f"helper of share {a.share}: F34 leaves {done} certified")


if __name__ == "__main__":
    main()
