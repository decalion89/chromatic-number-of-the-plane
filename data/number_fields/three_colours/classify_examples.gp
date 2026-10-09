read("classify.gp");  \\ run from this directory: gp -q classify_examples.gp
L = [["Q(sqrt2,sqrt7)", polcompositum(x^2-2, x^2-7)[1]], ["Q(sqrt2,sqrt31)", polcompositum(x^2-2, x^2-31)[1]], ["Q(sqrt10,sqrt38)", polcompositum(x^2-10, x^2-38)[1]], ["Q(sqrt2,sqrt55)", polcompositum(x^2-2, x^2-55)[1]], ["Q(sqrt3,sqrt5)", polcompositum(x^2-3, x^2-5)[1]], ["Q(sqrt3,sqrt7)", polcompositum(x^2-3, x^2-7)[1]], ["Q(sqrt2,sqrt3)", polcompositum(x^2-2, x^2-3)[1]], ["Q(c7,sqrt7)", polcompositum(x^3+x^2-2*x-1, x^2-7)[1]], ["Q(c7,sqrt2)", polcompositum(x^3+x^2-2*x-1, x^2-2)[1]], ["Q(sqrt11)", x^2-11], ["Q(sqrt7)", x^2-7], ["Q(sqrt3)", x^2-3], ["Q(sqrt2)", x^2-2], ["Q(cbrt2)", x^3-2]];
for(k = 1, #L, print(L[k][1], ": [a, b] = ", crit(L[k][2])));
quit;
