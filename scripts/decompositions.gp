\\ The decompositions of 2 and 3 used in notes/local_colourings.md (sections 8, 10 and 11) and in the
\\ note (papers/planes-4-chromatic), recomputed with PARI/GP.  Run from the root of the repository:
\\     gp -q scripts/decompositions.gp
\\ Each output line is: the generators of a multiquadratic field Q(sqrt a, sqrt b, ...), a prime p, and
\\ the sorted list of [e, f] (ramification index, residue degree) of the primes of the field over p.
\\ tests/test_decompositions.py runs this script and checks every line.

\\ A defining polynomial of Q(sqrt v[1], ..., sqrt v[k]); the square classes of v must be independent.
mq(v) = my(P = x^2 - v[1]); for (i = 2, #v, P = polcompositum(P, x^2 - v[i])[1]); P;

dec(v, p) = my(nf = nfinit(mq(v))); vecsort(apply(pr -> [pr.e, pr.f], idealprimedec(nf, p)));

show(v, p) = print(v, " ", p, " ", dec(v, p));

\\ The two theorems (sections 8 and 10): L and K = L(i).
show([3, 11], 2);
show([-1, 3, 11], 2);
show([2, 3], 2);
show([-1, 2, 3], 2);

\\ Q(sqrt3, sqrt q) for the primes q < 75 (section 11).
forprime (q = 2, 73, if (q != 3, show([3, q], 2); show([3, q], 3)));

\\ Nine fields with more square roots (section 11): the residue degree over 3 is 1 exactly when every
\\ generator other than 3 is 1 mod 3.
{
foreach ([[3, 7], [3, 7, 13], [3, 7, 13, 19], [3, 5], [3, 11], [3, 5, 7], [3, 7, 11], [3, 13, 17],
          [3, 7, 13, 23]], v, show(v, 3));
}
quit;
