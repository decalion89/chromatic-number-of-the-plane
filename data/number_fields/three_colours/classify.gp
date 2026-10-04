\\ classify.gp: Theorem B's criterion for F = Q[x]/(f): (a) a prime above 2 ramified in F(i); (b) a prime above 3 of
\\ residue degree 1.  Returns [a, b].  F is built in the variable y so that F(i) = F[x]/(x^2 + 1).
crit(f) = {
  my(K = nfinit(subst(f, x, y)), a = 0, b = 0, P3, P2, dd);
  if(#nfroots(K, x^2 + 1), return([0, 0]));          \\ i in F: neither condition can hold
  P3 = idealprimedec(K, 3); P2 = idealprimedec(K, 2);
  for(j = 1, #P3, if(P3[j].f == 1, b = 1));
  dd = rnfdisc(K, x^2 + 1)[1];
  for(j = 1, #P2, if(idealval(K, dd, P2[j]) > 0, a = 1));
  [a, b];
}
