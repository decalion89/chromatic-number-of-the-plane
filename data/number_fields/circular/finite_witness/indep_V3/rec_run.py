"""Run a script with pysat's Solver wrapped so that every add_clause / solve call is logged (sha256 of the stream).
usage: python3 rec_run.py LOG SCRIPT args...   (the script is run with runpy, sys.argv = [SCRIPT, args...])"""
import sys, hashlib, runpy
import pysat.solvers as ps

log = open(sys.argv[1], 'w')
_Orig = ps.Solver


class Rec(_Orig):
    _n = 0

    def __init__(self, *a, bootstrap_with=None, **k):
        Rec._n += 1
        self._id = Rec._n
        log.write(f'new solver {self._id} {k.get("name", a[0] if a else "")}\n')
        super().__init__(*a, **k)
        if bootstrap_with:
            for c in bootstrap_with:
                self.add_clause(c)

    def add_clause(self, c, no_return=True):
        log.write(f'{self._id} c ' + ' '.join(map(str, c)) + '\n')
        return super().add_clause(c, no_return)

    def solve(self, assumptions=[]):
        r = super().solve(assumptions=assumptions)
        log.write(f'{self._id} solve {len(assumptions)} -> {r}\n')
        return r

    def solve_limited(self, assumptions=[], expect_interrupt=False):
        r = super().solve_limited(assumptions=assumptions, expect_interrupt=expect_interrupt)
        log.write(f'{self._id} solve_limited {len(assumptions)} -> {r}\n')
        return r


ps.Solver = Rec
script = sys.argv[2]
sys.argv = sys.argv[2:]
try:
    runpy.run_path(script, run_name='__main__')
finally:
    log.close()
