\\ Theorem A of notes/local_global.md: for a number field F, chi(F^2) = 2 iff some prime of F above 2 ramifies in
\\ F(i). ram2(f) returns 1 when that holds for the field F = Q[y]/(f(y)), using the relative discriminant of F(i)/F.
\\ Run: gp -q two_colour_criterion.gp < /dev/null
ram2(f) = {
  my(F = nfinit(f), dd, PF);
  if (#nfroots(F, x^2 + 1), return(0));          \\ i in F: F(i) = F and nothing ramifies
  dd = rnfdisc(F, x^2 + 1)[1];                    \\ relative discriminant ideal of F(i)/F
  PF = idealprimedec(F, 2);
  for (j = 1, #PF, if (idealval(F, dd, PF[j]) > 0, return(1)));
  0
}
quad(d) = ram2(y^2 - d);
biq(a, b) = ram2(polcompositum(y^2 - a, y^2 - b)[1]);
{
  my(bad = 0, n = 0);
  forsquarefree(d = 2, 400, my(m = d[1]); n++; if (quad(m) != (m % 4 != 3), bad++; print("MISMATCH d=", m)));
  print("quadratic: ", n, " fields d <= 400, mismatches with chi = 2 iff d != 3 mod 4 (Johnson, Fischer): ", bad);
}
{
  my(bad = 0, L = [y^3 - 2, y^3 - 3, y^3 - y - 1, y^3 - 3*y - 1, y^3 + y^2 - 2*y - 1, y^5 - 2, y^3 - 5, y^5 - y - 1, y^7 - 3]);
  for (j = 1, #L, if (ram2(L[j]) != 1, bad++; print("MISMATCH odd degree ", L[j])));
  print("odd degree: ", #L, " fields, mismatches with chi = 2 (Moorhouse, Theorem 7.1): ", bad);
}
{
  my(n = 0, bad = 0);
  for (a = 2, 60, if (!issquarefree(a), next); for (b = a + 1, 120, if (!issquarefree(b) || gcd(a, b) != 1, next);
     n++; my(three = (a % 4 == 3) || (b % 4 == 3) || (core(a * b) % 4 == 3));
     if (biq(a, b) != !three, bad++; print("MISMATCH biquadratic ", [a, b]))));
  print("biquadratic: ", n, " fields Q(sqrt a, sqrt b), a <= 60, a < b <= 120 coprime; mismatches with chi = 2 iff no quadratic subfield has c = 3 mod 4: ", bad);
}
