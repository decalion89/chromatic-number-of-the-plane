# Test: for Cayley graphs Cay(Z/n x Z/m, S), S = -S, is  chi <= 3  <=>  exists character xi with xi(S) in [1/3, 2/3] (mod 1)?
import itertools, random, sys
from pysat.solvers import Solver
def three_col(n, m, S):
    N = n*m; idx = lambda a, b: (a % n)*m + (b % m)
    v = lambda x, c: 3*x + c + 1
    s = Solver(name="cadical195")
    for x in range(N): s.add_clause([v(x, c) for c in range(3)])
    for a in range(n):
        for b in range(m):
            for (p, q) in S:
                y = idx(a+p, b+q); x = idx(a, b)
                if x < y:
                    for c in range(3): s.add_clause([-v(x, c), -v(y, c)])
                elif x == y: return False
    return s.solve()
def char_ok(n, m, S):
    for j in range(n):
        for k in range(m):
            ok = True
            for (p, q) in S:
                t = (p*j*m + q*k*n) % (n*m)          # xi(p,q) = p j/n + q k/m, numerator over n m
                if not (n*m <= 3*t <= 2*n*m): ok = False; break
            if ok: return True
    return False
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
bad = 0; cnt = {}
for trial in range(4000):
    n = random.randint(3, 13); m = random.choice([1, 1, 2, 3, 4, 5, 6])
    r = random.randint(1, 3)
    S = set()
    for _ in range(r):
        p, q = random.randrange(n), random.randrange(m)
        if (p, q) == (0, 0): continue
        S.add((p, q)); S.add(((-p) % n, (-q) % m))
    if not S: continue
    S = sorted(S)
    a = three_col(n, m, S); b = char_ok(n, m, S)
    cnt[(a, b)] = cnt.get((a, b), 0) + 1
    if a != b:
        bad += 1
        if bad < 10: print("MISMATCH", n, m, S, "3col", a, "char", b)
print(cnt, "mismatches", bad)
