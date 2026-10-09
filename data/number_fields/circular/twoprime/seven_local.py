"""Local facts at 7 (exact, F_49 = Z[i]/7):
 E7 = { e in F_49 : Re(conj(e) z) in {2,3,4,5} for every z in mu_8 = {z : z conj(z) = 1} }  (e = 7*epsilon);
 checks |E7| = 8, E7 is one mu_8-orbit, 0 not in E7, E7 contains no coset of a non-zero additive subgroup of F_49,
 the first coordinate maps E7 into {2,3,4,5}, and lists E7 and the 7-torsion points c = 65*eps' of the two-prime
 window (eps' = 65^{-1} eps mod 7 Z[i])."""
p = 7
F = [(a, b) for a in range(p) for b in range(p)]
def mul(x, y): return ((x[0]*y[0] - x[1]*y[1]) % p, (x[0]*y[1] + x[1]*y[0]) % p)
def conj(x): return (x[0], (-x[1]) % p)
mu8 = [z for z in F if mul(z, conj(z)) == (1, 0)]
assert len(mu8) == 8
E7 = [e for e in F if all(mul(conj(e), z)[0] in (2, 3, 4, 5) for z in mu8)]
print("mu_8 =", mu8)
print("E7 (7*epsilon mod 7) =", E7, " size", len(E7))
orbit = sorted({mul(E7[0], z) for z in mu8})
print("one mu_8-orbit:", orbit == sorted(E7))
print("0 in E7:", (0, 0) in E7, "; first coordinates:", sorted({e[0] for e in E7}))
# cosets of lines: c + t d, t in F_7, d != 0
bad = []
for d in F:
    if d == (0, 0): continue
    for c in F:
        line = {((c[0] + t * d[0]) % p, (c[1] + t * d[1]) % p) for t in range(p)}
        if line <= set(E7): bad.append((c, d))
print("cosets of non-zero subgroups inside E7:", bad if bad else "none")
inv65 = pow(65, -1, p)
pts = sorted(((e[0] * inv65) % p, (e[1] * inv65) % p) for e in E7)
print("classes c/65 of the 7-torsion points of S(1,1) (expected (1,2),(1,5),(2,1),(2,6),(5,1),(5,6),(6,2),(6,5)):", pts)
s = (sum(e[0] for e in E7) % p, sum(e[1] for e in E7) % p)
print("sum of E7:", s)
