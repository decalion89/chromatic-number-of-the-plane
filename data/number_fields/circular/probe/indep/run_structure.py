"""Check Proposition D1 by exact computation of S_N^(r), N = 5^k, k <= K, for given r values.
usage: python3 run_structure.py K r1 r2 ...   (r as fractions like 17/56)"""
import sys
import time
from fractions import Fraction as Fr
from sn import lift_levels, direct_level, classify, make_refs, GN_brute, GN_rho, kappa, in_S, cstar, QN, Pk, Yk


def canon_mod(p, N):
    return (p[0] % N, p[1] % N)


def piece_key(P, N):
    # canonical description modulo N: sorted vertices shifted so that the min vertex lies in [0,N)^2
    V = sorted(P.V)
    v0 = V[0]
    t = (v0[0] - (v0[0] % N), v0[1] - (v0[1] % N))
    return tuple((v[0] - t[0], v[1] - t[1]) for v in V)


def main():
    K = int(sys.argv[1])
    rs = [Fr(a) for a in sys.argv[2:]]
    direct_upto = int(__import__("os").environ.get("DIRECT", "2"))
    brute_upto = int(__import__("os").environ.get("BRUTE", "3"))
    for r in rs:
        s = Fr(1, 2) - r
        eps = Fr(1, 3) - r
        t0 = time.time()
        levels = lift_levels(K, r)
        print(f"=== r = {r} (= {float(r):.7f}), s = {s}, eps = {eps}; lifting took {time.time() - t0:.1f}s")
        for k in range(1, K + 1):
            N = 5 ** k
            pieces = levels[k]
            refs = make_refs(k, r)
            labs = [classify(P, k, r, refs) for P in pieces]
            nc = sum(1 for l in labs if l and l[0] == "c")
            nq = sorted(str(l[1]) for l in labs if l and l[0] == "q")
            junk = [P for P, l in zip(pieces, labs) if l is None]
            ndeg = sum(1 for P in junk if len(P.V) == 1)
            line = f"  k={k} N={N}: {len(pieces)} pieces; c-component x{nc}; q-components x{len(nq)} (distinct q: {len(set(nq))}); junk {len(junk)} ({ndeg} isolated points)"
            print(line)
            # radius checks on the reference polygons
            Pv = refs[0][2]
            maxP = max(v[0] ** 2 + v[1] ** 2 for v in Pv)
            ok_P = maxP <= Fr(10, 9) * s * s
            maxY = max(max(v[0] ** 2 + v[1] ** 2 for v in ref) for lab, T, ref in refs[1:])
            ok_Y = maxY <= 2 * eps * eps
            print(f"     max|x|^2 on P_k = {maxP} <= 10 s^2/9 = {Fr(10, 9) * s * s}: {ok_P};  max|y|^2 on Y_k = {maxY} <= 2 eps^2 = {2 * eps * eps}: {ok_Y};  #vertices P_k: {len(Pv)}")
            # pointwise sanity: every vertex and vertex-centroid of every piece lies in S_N^(r) (brute-force G_N)
            if k <= brute_upto:
                G = GN_brute(N)
                assert len(G) == 4 * (2 * k + 1), len(G)
                assert sorted(G) == sorted(GN_rho(k))
                bad = 0
                for P in pieces:
                    for v in list(P.V) + [P.centroid_v()]:
                        if not in_S(v, G, r):
                            bad += 1
                print(f"     brute-force G_N ({len(G)} rotations, equal to rho^j form): vertices+centroids not in S: {bad}")
            if junk:
                for P in junk[:40]:
                    G = GN_brute(N) if k <= brute_upto else GN_rho(k)
                    v = P.centroid_v()
                    print(f"     junk: {len(P.V)} vertices, centroid ({v[0]}, {v[1]}), kappa(centroid) = {kappa(v, G)}")
            if k <= direct_upto:
                t1 = time.time()
                D = direct_level(k, r)
                A = sorted(piece_key(P, N) for P in pieces)
                B = sorted(piece_key(P, N) for P in D)
                print(f"     direct enumeration ({time.time() - t1:.1f}s): {len(D)} pieces; identical to lifting: {A == B}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
