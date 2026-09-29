\\ gpcheck.gp -- independent exact check (PARI/gp) of a graph given by V (rows: 7*coefficients of 1,z,...,z^11,
\\ z = zeta21) and E (1-based index pairs). Checks: vertices distinct; every edge has d*conj(d) = 1 in Q(zeta21).
read("g1023_gp_input.gp");
f = polcyclo(21);
el(r) = Mod(sum(k = 1, 12, r[k] * x^(k-1)) / 7, f);
cj(a) = Mod(subst(lift(a), x, x^20), f);
n = matsize(V)[1]; m = matsize(E)[1];
X = vector(n, i, el(V[i,]));
distinct = (#Set(vector(n, i, lift(X[i]))) == n);
ok = 1; for (k = 1, m, my(d = X[E[k,2]] - X[E[k,1]]); if (d * cj(d) != 1, ok = 0; print("bad edge ", k)));
print("PARI/gp: n = ", n, ", m = ", m, "; vertices distinct: ", distinct, "; all edges exactly unit (d*conj(d) = 1): ", ok);
\\ numeric embedding z -> exp(2 pi i/21): count all unit-distance pairs among the vertices (float, 1e-20)
\p 60
zz = exp(2*Pi*I/21);
Z = vector(n, i, subst(lift(X[i]), x, zz));
cnt = 0; for (i = 1, n, for (j = i+1, n, if (abs(abs(Z[i]-Z[j]) - 1) < 1e-20, cnt++)));
print("PARI/gp: unit-distance pairs among the vertices (numeric, 60 digits): ", cnt);
quit;
