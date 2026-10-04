\\ Independent check (PARI/GP only) of the N -> infinity limit test for a finite set V of unit vectors of F^2,
\\ F = Q[a,b]/(fa(a), fb(b)).  phi(v) = sum_j c_j X_j(v) (X_j(v) in Q(i): coefficient of the j-th monomial in X + iY).
\\ Feasible iff some z in L = (column space of the real matrix) cap Z^{2m} has, for every v, 6phi(v) = z_v with
\\ z_v = (1,1) mod 2 and (0,0) mod 3 [type c] or z_v = (0,0) mod 2 and both coords nonzero mod 3 [type q].
red(x) = lift(lift(Mod(Mod(x, fb), fa)));
coef(x, i, j) = polcoef(polcoef(red(x), i, a), j, b);
m = #V;
for(k = 1, m, if(red(V[k][1]^2 + V[k][2]^2 - 1) != 0, error("not a unit vector: ", k)));
nb = da * db;
A = matrix(2*m, 2*nb);
{
for(k = 1, m,
  my(X = V[k][1], Y = V[k][2], col = 0);
  for(i = 0, da - 1, for(j = 0, db - 1,
    my(x = coef(X, i, j), y = coef(Y, i, j));
    col++;
    \\ unknowns (cr, ci) for monomial col: Re(c*(x+iy)) = cr*x - ci*y ; Im = cr*y + ci*x
    A[2*k-1, 2*col-1] = x; A[2*k-1, 2*col] = -y;
    A[2*k, 2*col-1] = y; A[2*k, 2*col] = x)));
}
D = denominator(A); Ai = A * D;
K = matkerint(Ai~);            \\ integer basis of {y : y~ * Ai = 0}
L = if(#K, matkerint(K~), matid(2*m));   \\ saturated lattice {z : K~ z = 0}
r = #L;
print("m = ", m, ", rank of L = ", r, " (expected ", 2*nb, ")");
if(r != 2*nb, error("unexpected rank"));
\\ mod 2: all patterns
pats = List();
{
forvec(e = vector(r, i, [0, 1]),
  my(x = (L * e~) % 2, ok = 1, p = vector(m));
  for(k = 1, m,
    my(s = x[2*k-1], t = x[2*k]);
    if(s == 1 && t == 1, p[k] = 1, if(s == 0 && t == 0, p[k] = 0, ok = 0; break)));
  if(ok, listput(pats, p)));
}
pats = Set(Vec(pats));
print("surviving mod-2 patterns: ", #pats);
feas = 0;
{
forvec(e = vector(r, i, [0, 2]),
  my(y = (L * e~) % 3);
  for(n = 1, #pats,
    my(p = pats[n], ok = 1);
    for(k = 1, m,
      my(s = y[2*k-1], t = y[2*k]);
      if(p[k] == 1, if(s != 0 || t != 0, ok = 0; break), if(s == 0 || t == 0, ok = 0; break)));
    if(ok, feas = 1; print("FEASIBLE, pattern ", p); break));
  if(feas, break));
}
if(!feas, print("INFEASIBLE"));
quit;
