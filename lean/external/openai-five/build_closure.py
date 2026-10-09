"""Compile the import closure of one module of the openai/math Lean library against an existing Mathlib build.

usage: python3 build_closure.py OAI_LEAN_DIR LEAN_BIN LEAN_PATH OUT_DIR [JOBS] [TARGET]

OAI_LEAN_DIR  the lean/ directory of a checkout of https://github.com/openai/math (it holds OAI/...)
LEAN_BIN      the lean binary of toolchain v4.34.1
LEAN_PATH     colon-separated .lake/build/lib/lean directories of a Mathlib build at commit d13f23b7
              (mathlib, batteries, aesop, Qq, plausible, proofwidgets, LeanSearchClient, importGraph),
              for example that of lean/ in this repository after `lake build`
OUT_DIR       where the .olean files go
TARGET        default OAI.Geometry.PlaneColoring.Five

Modules are compiled in dependency order, JOBS at a time (default 2). One line per module is printed:
`ok M s`, `FAILED M`, or `SKIPPED M`; the compiler's messages go to OUT_DIR/M.out (all empty when nothing is
reported). Then compile PlaneSix.lean with LEAN_PATH=OUT_DIR:LEAN_PATH to see the restated theorem and its axioms.
"""
import os, re, subprocess, sys, time

src, lean, lp, out = sys.argv[1:5]
jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 2
target = sys.argv[6] if len(sys.argv) > 6 else 'OAI.Geometry.PlaneColoring.Five'
path = lambda m: os.path.join(*m.split('.')) + '.lean'
olean = lambda m: os.path.join(out, *m.split('.')) + '.olean'
deps = {}
def visit(m):
    if m in deps:
        return
    deps[m] = []
    for line in open(os.path.join(src, path(m)), encoding='utf-8'):
        mm = re.match(r'\s*(public\s+)?import\s+(\S+)', line)
        if mm and mm.group(2).startswith('OAI'):
            deps[m].append(mm.group(2))
    for d in deps[m]:
        visit(d)
visit(target)
print(f'{len(deps)} modules, {sum(sum(1 for _ in open(os.path.join(src, path(m)), encoding="utf-8")) for m in deps)} lines',
      flush=True)
done, failed, running = set(), set(), {}
t0 = time.time()
while len(done) + len(failed) < len(deps):
    for m, (p, ts) in list(running.items()):
        if p.poll() is not None:
            del running[m]
            if p.returncode == 0 and os.path.exists(olean(m)):
                done.add(m); print(f'ok {m} {time.time() - ts:.0f}s', flush=True)
            else:
                failed.add(m); print(f'FAILED {m} rc={p.returncode}', flush=True)
    for m in deps:
        if m not in done | failed and m not in running and any(d in failed for d in deps[m]):
            failed.add(m); print(f'SKIPPED {m}', flush=True)
    ready = [m for m in deps if m not in done | failed and m not in running and all(d in done for d in deps[m])]
    while ready and len(running) < jobs:
        m = ready.pop(0)
        os.makedirs(os.path.dirname(olean(m)), exist_ok=True)
        log = open(olean(m)[:-6] + '.out', 'w')
        p = subprocess.Popen([lean, '-o', olean(m), '-i', olean(m)[:-6] + '.ilean', path(m)], cwd=src,
                             env=dict(os.environ, LEAN_PATH=out + ':' + lp), stdout=log, stderr=subprocess.STDOUT)
        running[m] = (p, time.time())
    if not running and not ready:
        break
    time.sleep(2)
print(f'END built {len(done)}, failed {len(failed)}, {time.time() - t0:.0f}s', flush=True)
sys.exit(1 if failed else 0)
