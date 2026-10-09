/* Referee program (method: exact vertex enumeration of the line arrangement).
 *
 * S^r(N) = { c = x+iy : Re(conj(c) g) in [r, 1-r] + Z for every rotation g with N g in Z[i] }.
 * For N = 5^k 13^m these rotations are exactly G(k,m).  S^r(N) is invariant under N Z[i].
 * Every rotation g = (A+Bi)/N gives the functionals Re(conj(c) g) = (A x + B y)/N and
 * Im(conj(c) g) = (B x - A y)/N; up to sign they are listed once.
 *
 * Every component of S^r(N) (a fixed vector of strip indices) is a non-empty bounded convex polygon, so it has a
 * vertex, which is the intersection of two boundary lines  f_k = n + r  or  f_k = n + 1 - r  of two non-parallel
 * functionals, and which lies in S^r(N).  Every component lies inside a square [a+r, a+1-r] x [b+r, b+1-r], so
 * modulo N Z[i] it has a translate inside [0,N)^2.  This program enumerates ALL intersection points of pairs of
 * boundary lines in [0,N)^2 and prints those that lie in S^r(N) (exact integer arithmetic, 128-bit).
 *
 * r = rn/rd.  A point is printed as  Xn Yn D  meaning  x = N*Xn/D, y = N*Yn/D  (D > 0).
 * usage: vertex_enum N rn rd
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef long long ll;

static ll FA[4096], FB[4096];
static int nf = 0;

static i128 fmod128(i128 a, i128 m) { /* m > 0, result in [0, m) */
    i128 q = a % m;
    if (q < 0) q += m;
    return q;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: vertex_enum N rn rd\n"); return 1; }
    ll N = atoll(argv[1]), rn = atoll(argv[2]), rd = atoll(argv[3]);
    if (!(0 < 2 * rn && 2 * rn < rd)) { fprintf(stderr, "need 0 < r < 1/2\n"); return 1; }
    /* rotations */
    int nrot = 0;
    for (ll A = -N; A <= N; A++) {
        ll B2 = N * N - A * A;
        ll B = 0;
        while ((B + 1) * (B + 1) <= B2) B++;
        if (B * B != B2) continue;
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            ll Bs = sgn * B;
            if (B == 0 && sgn == 1) continue; /* avoid duplicate when B = 0 */
            nrot++;
            ll cand[2][2] = {{A, Bs}, {Bs, -A}};
            for (int t = 0; t < 2; t++) {
                ll a = cand[t][0], b = cand[t][1];
                if (a < 0 || (a == 0 && b < 0)) { a = -a; b = -b; }
                int dup = 0;
                for (int u = 0; u < nf; u++) if (FA[u] == a && FB[u] == b) { dup = 1; break; }
                if (!dup) { FA[nf] = a; FB[nf] = b; nf++; }
            }
        }
    }
    /* put the functionals of g = 1 first: (N,0) and (0,N) */
    for (int u = 0; u < nf; u++) {
        if (FA[u] == N && FB[u] == 0) { ll ta = FA[0], tb = FB[0]; FA[0] = FA[u]; FB[0] = FB[u]; FA[u] = ta; FB[u] = tb; }
    }
    for (int u = 0; u < nf; u++) {
        if (FA[u] == 0 && FB[u] == N) { ll ta = FA[1], tb = FB[1]; FA[1] = FA[u]; FB[1] = FB[u]; FA[u] = ta; FB[u] = tb; }
    }
    fprintf(stderr, "N = %lld, r = %lld/%lld: %d rotations, %d functionals up to sign\n", N, rn, rd, nrot, nf);
    printf("# N %lld r %lld %lld rotations %d functionals %d\n", N, rn, rd, nrot, nf);
    for (int u = 0; u < nf; u++) printf("# F %lld %lld\n", FA[u], FB[u]);
    /* boundary levels of each functional on [0,N]^2: value v = (A x + B y)/N in [lo, hi] */
    ll *lev[4096];
    int nlev[4096];
    for (int u = 0; u < nf; u++) {
        ll lo = (FA[u] < 0 ? FA[u] : 0) + (FB[u] < 0 ? FB[u] : 0);
        ll hi = (FA[u] > 0 ? FA[u] : 0) + (FB[u] > 0 ? FB[u] : 0);
        int cnt = 0;
        lev[u] = (ll *)malloc(sizeof(ll) * (size_t)(2 * (hi - lo + 4)));
        for (ll n = lo - 1; n <= hi + 1; n++) {
            lev[u][cnt++] = rd * n + rn;          /* f = n + r      (level L: f = L/rd) */
            lev[u][cnt++] = rd * n + rd - rn;     /* f = n + 1 - r */
        }
        nlev[u] = cnt;
    }
    long long npairs = 0, ninside = 0, nout = 0;
    for (int k = 0; k < nf; k++) {
        for (int l = k + 1; l < nf; l++) {
            ll det = FA[k] * FB[l] - FA[l] * FB[k];
            if (det == 0) { fprintf(stderr, "parallel functionals?\n"); return 2; }
            i128 D0 = (i128)rd * det;
            int sg = D0 < 0 ? -1 : 1;
            i128 D = D0 * sg;
            ll adet = det < 0 ? -det : det;
            i128 lowq = (i128)rn * adet, highq = (i128)(rd - rn) * adet;
            for (int a = 0; a < nlev[k]; a++) {
                ll Lk = lev[k][a];
                for (int b = 0; b < nlev[l]; b++) {
                    ll Ll = lev[l][b];
                    npairs++;
                    i128 Xn = ((i128)Lk * FB[l] - (i128)Ll * FB[k]) * sg;
                    if (Xn < 0 || Xn >= D) continue;
                    i128 Yn = ((i128)FA[k] * Ll - (i128)FA[l] * Lk) * sg;
                    if (Yn < 0 || Yn >= D) continue;
                    ninside++;
                    int ok = 1;
                    for (int m = 0; m < nf; m++) {
                        i128 Q = (i128)FA[m] * Xn + (i128)FB[m] * Yn;
                        i128 q = fmod128(Q, D);
                        if (q < lowq || q > highq) { ok = 0; break; }
                    }
                    if (ok) {
                        nout++;
                        printf("%lld %lld %lld\n", (ll)Xn, (ll)Yn, (ll)D);
                    }
                }
            }
        }
    }
    fprintf(stderr, "line pairs tested %lld, intersection points in [0,N)^2: %lld, in S^r: %lld\n", npairs, ninside, nout);
    printf("# pairs %lld inside %lld out %lld\n", npairs, ninside, nout);
    return 0;
}
