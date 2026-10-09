\\ The new half of Theorem A (notes/local_global.md) on fields where no quadratic subfield explains it.
\\ For F = Q[y]/(f), unit vectors of F^2 are u = z/zbar = ((a^2 - b^2)/N, 2ab/N), N = a^2 + b^2, for z = a + ib in
\\ O_{F(i)} (Hilbert 90). An integer relation sum c_k u_k = 0 with sum c_k odd is a closed walk of odd length
\\ sum |c_k|, so chi(F^2) >= 3. Fields where some prime above 2 ramifies in F(i) (chi = 2) can only give even sums.
\\ Run: gp -q odd_walks.gp < /dev/null
ram2(f) = {
  my(F = nfinit(f), dd, PF);
  if (#nfroots(F, x^2 + 1), return(0));
  dd = rnfdisc(F, x^2 + 1)[1];
  PF = idealprimedec(F, 2);
  for (j = 1, #PF, if (idealval(F, dd, PF[j]) > 0, return(1)));
  0
}
\\ an integral basis of O_L, L = F[x]/(x^2 + 1), as pairs [a, b] (a + i b) with a, b in F
intbasis_rel(f) = {
  my(F = nfinit(f), rnf = rnfinit(F, x^2 + 1), Labs = nfinit(rnf.polabs), B = Labs.zk, out = List());
  for (k = 1, #B,
    my(r = lift(rnfeltabstorel(rnf, Mod(B[k], rnf.polabs))));
    listput(out, [polcoef(r, 0, x), polcoef(r, 1, x)]));
  [F, Vec(out)]
}
\\ m random unit vectors (coefficients in [-B, B] on the integral basis, seed s); returns
\\ [rank of the relation lattice, number of basis relations with odd sum, shortest odd walk found]
oddwalks(f, m, B, s) = {
  my(tmp = intbasis_rel(f), F = tmp[1], IB = tmp[2], U = List(), cols = List(), M, K, nodd = 0, best = 0);
  setrand(s);
  while (#U < m,
    my(a = 0, b = 0);
    for (k = 1, #IB, my(c = random(2 * B + 1) - B); a += c * IB[k][1]; b += c * IB[k][2]);
    a = Mod(lift(a), f); b = Mod(lift(b), f);
    my(N = a^2 + b^2);
    if (N == 0, next);
    my(u = [(a^2 - b^2) / N, 2 * a * b / N]);
    if (u[1]^2 + u[2]^2 != 1, error("not a unit vector"));
    listput(U, u);
    listput(cols, concat(nfalgtobasis(F, lift(u[1])), nfalgtobasis(F, lift(u[2])))));
  M = matconcat(Vec(cols)); M = M * denominator(M);
  K = matkerint(M);
  for (j = 1, #K,
    my(c = K[, j]);
    if (sum(k = 1, m, c[k] * U[k][1]) != 0 || sum(k = 1, m, c[k] * U[k][2]) != 0, error("relation check failed"));
    if (vecsum(c) % 2,
      nodd++;
      my(len = sum(k = 1, m, abs(c[k])));
      if (!best || len < best, best = len)));
  [#K, nodd, best]
}
{
  my(odd = [y^4 - 2*y^3 - 2*y - 3, y^4 - 6*y^2 + 4*y + 2, y^4 - 5*y^2 - 7],
     even = [polcompositum(y^2 - 2, y^2 - 5)[1], y^4 - 2, y^3 - 2]);
  for (k = 1, #odd, my(f = odd[k], r = oddwalks(f, 40, 1, 7));
    print("odd walk: ", f, ", quadratic subfields ", #nfsubfields(f, 2), ", ram2 ", ram2(f),
          ", relations ", r[1], ", odd ", r[2], ", shortest odd walk found ", r[3]));
  for (k = 1, #even, my(f = even[k], r = oddwalks(f, 40, 1, 7));
    print("even only: ", f, ", ram2 ", ram2(f), ", relations ", r[1], ", odd ", r[2]));
}
