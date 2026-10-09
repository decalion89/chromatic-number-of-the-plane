"""For the 140 vectors U of Theorem C: every 7-adic character z -> Re(conj(e) * (z mod w)) / 7, e in A'_7,
at every place w above 7, keeps U in {2,3,4,5}/7, so it is an optimal character if kappa(U) = 2/7.
Lemma 22's key step then predicts that its tight set T = {u : value = 2/7} carries a positive integer relation
(0 in the convex hull of T in L (x) R = R^4).  If it did not, a perturbation would give kappa(U) > 2/7,
contradicting the certificate of Theorem C.  We find such relations and verify them exactly."""
from fractions import Fraction as Fr
import numpy as np, sympy, json
from scipy.optimize import linprog
from thmC_vectors import build_U

def red(fr): return fr.numerator*pow(fr.denominator,-1,7)%7
A7=[(2,3),(2,4),(3,2),(3,5),(4,2),(4,5),(5,3),(5,4)]
out={}
for d in (11,35):
    U=build_U(d)
    roots=[2,5] if d==11 else [0]
    res=[]
    for root in roots:
        for (ea,eb) in A7:
            vals=[]
            for (a,b,c,e) in U:
                x=(red(a)+root*red(b))%7; y=(red(c)+root*red(e))%7
                # Re(conj(e) z) with e = ea + eb i, z = x + y i : ea*x + eb*y
                vals.append((ea*x+eb*y)%7)
            assert set(vals) <= {2,3,4,5}, (d,root,ea,eb,set(vals))
            T=[U[i] for i,v in enumerate(vals) if v==2]
            # LP: lambda >= 0, sum = 1, sum lambda_u u = 0
            Am=np.array([[float(q) for q in u] for u in T]).T
            Aeq=np.vstack([Am,np.ones(len(T))]); beq=np.array([0,0,0,0,1.0])
            lp=linprog(np.zeros(len(T)),A_eq=Aeq,b_eq=beq,bounds=[(0,None)]*len(T),method='highs')
            ok=False; supp=None
            if lp.status==0:
                supp=[i for i,l in enumerate(lp.x) if l>1e-9]
                # exact: find positive rational solution on the support
                M=sympy.Matrix([[u[k] for u in [T[i] for i in supp]] for k in range(4)])
                ns=M.nullspace()
                # try to find a strictly positive combination: use LP solution direction rounded
                lam=[sympy.Rational(lp.x[i]).limit_denominator(10**6) for i in supp]
                # project to nullspace exactly
                if ns:
                    B=sympy.Matrix.hstack(*ns)
                    coef=(B.T*B).inv()*B.T*sympy.Matrix(lam)
                    lamx=B*coef
                    if all(v>0 for v in lamx):
                        den=sympy.ilcm(*[sympy.fraction(v)[1] for v in lamx])
                        ints=[int(v*den) for v in lamx]
                        g=sympy.igcd(*ints); ints=[v//g for v in ints]
                        tot=[sum(n*T[i][k] for n,i in zip(ints,supp)) for k in range(4)]
                        ok=all(t==0 for t in tot) and all(n>0 for n in ints)
                        supp=(len(supp),max(ints),sum(ints))
            res.append(dict(place_root=root,e=(ea,eb),tight=len(T),relation_found=ok,support_maxcoef_length=supp))
    out[d]=res
    print(d, 'characters:',len(res),' all values in {2..5}/7: True;  tight-set sizes:',sorted(set(r['tight'] for r in res)),
          ' positive relation found & verified exactly for all:', all(r['relation_found'] for r in res), flush=True)
    print('   examples:', res[:3])
json.dump(out,open('thmC_tight7_results.json','w'),indent=1,default=str)
