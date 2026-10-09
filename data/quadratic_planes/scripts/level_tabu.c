/* level_tabu.c -- tabu search for a proper K-colouring of the level-r plane Cay((Z/p^r)^2, U_r), U_r = {u : N(u) = 1},
   N(a,b) = a^2 + b^2, for p = 3 (mod 4).  Neighbours are computed from the unit vectors, edges are not stored.
   usage: level_tabu p r K maxsec seed [init.txt]   (build: gcc -O3 -o level_tabu level_tabu.c) ; prints best conflict count over time and writes best colouring.
   Move rule (Hertz-de Werra tabucol with a sampled candidate list): pick the best (vertex, colour) move among
   all conflicting vertices if few, else among a random sample; tabu tenure L + 0.6*conflicts. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static unsigned long long rs = 88172645463325252ULL;
static inline unsigned long long xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: level_tabu p r K maxsec seed [init]\n"); return 1; }
    int p = atoi(argv[1]), r = atoi(argv[2]), K = atoi(argv[3]); double maxsec = atof(argv[4]);
    rs ^= strtoull(argv[5], 0, 10) * 2654435761ULL + 1;
    long M = 1; for (int i = 0; i < r; i++) M *= p;
    long n = M * M;
    /* unit vectors */
    int nu = 0; int *ua = malloc(sizeof(int) * 4 * M), *ub = malloc(sizeof(int) * 4 * M);
    for (long a = 0; a < M; a++) for (long b = 0; b < M; b++) if ((a * a + b * b) % M == 1 % M) { ua[nu] = a; ub[nu] = b; nu++; }
    fprintf(stderr, "p=%d r=%d: %ld vertices, %d unit vectors, K=%d\n", p, r, n, nu, K);
    unsigned char *col = malloc(n);
    int *g = calloc((size_t)n * K, sizeof(int));
    if (argc > 6) { FILE *f = fopen(argv[6], "r"); for (long v = 0; v < n; v++) { int c; if (fscanf(f, "%d", &c) != 1) return 2; col[v] = c % K; } fclose(f); }
    else for (long v = 0; v < n; v++) col[v] = xr() % K;
    for (long v = 0; v < n; v++) { long a = v / M, b = v % M;
        for (int k = 0; k < nu; k++) { long w = ((a + ua[k]) % M) * M + (b + ub[k]) % M; g[w * K + col[v]]++; } }
    long conf = 0; for (long v = 0; v < n; v++) conf += g[v * K + col[v]]; conf /= 2;
    /* list of conflicting vertices */
    long *lst = malloc(sizeof(long) * n), *pos = malloc(sizeof(long) * n), nl = 0;
    for (long v = 0; v < n; v++) { pos[v] = -1; if (g[v * K + col[v]] > 0) { pos[v] = nl; lst[nl++] = v; } }
    long *tabu = calloc((size_t)n * K, sizeof(long));
    long best = conf; long it = 0; time_t t0 = time(0), tl = t0;
    unsigned char *bestcol = malloc(n); memcpy(bestcol, col, n);
    fprintf(stderr, "start conflicts %ld\n", conf);
    while (conf > 0) {
        it++;
        /* candidate moves: sample up to S conflicting vertices */
        int S = nl < 400 ? (int)nl : 400;
        long bv = -1; int bc = -1, bd = 1 << 30, cnt = 0;
        for (int s = 0; s < S; s++) {
            long v = nl < 400 ? lst[s] : lst[xr() % nl];
            int cv = col[v], gv = g[v * K + cv];
            for (int c = 0; c < K; c++) if (c != cv) {
                int d = g[v * K + c] - gv;
                int istabu = tabu[v * K + c] > it;
                if (istabu && !(conf + d < best)) continue;   /* aspiration */
                if (d < bd) { bd = d; bv = v; bc = c; cnt = 1; }
                else if (d == bd && (xr() % (++cnt)) == 0) { bv = v; bc = c; }
            }
        }
        if (bv < 0) continue;
        int oc = col[bv]; col[bv] = bc; conf += bd;
        tabu[bv * K + oc] = it + 10 + (long)(0.6 * nl) + (xr() % 10);
        long a = bv / M, b = bv % M;
        for (int k = 0; k < nu; k++) {
            long w = ((a + ua[k]) % M) * M + (b + ub[k]) % M;
            g[w * K + oc]--; g[w * K + bc]++;
            int cw = col[w];
            int was = (cw == oc && g[w * K + oc] + 1 > 0), now = g[w * K + cw] > 0;
            (void)was;
            if (now && pos[w] < 0) { pos[w] = nl; lst[nl++] = w; }
            else if (!now && pos[w] >= 0) { long q = pos[w]; lst[q] = lst[--nl]; pos[lst[q]] = q; pos[w] = -1; }
        }
        int nowv = g[bv * K + bc] > 0;
        if (nowv && pos[bv] < 0) { pos[bv] = nl; lst[nl++] = bv; }
        else if (!nowv && pos[bv] >= 0) { long q = pos[bv]; lst[q] = lst[--nl]; pos[lst[q]] = q; pos[bv] = -1; }
        if (conf < best) { best = conf; memcpy(bestcol, col, n); }
        if ((it & 1023) == 0) {
            time_t t = time(0);
            if (t - tl >= 60) { fprintf(stderr, "t=%lds it=%ld conflicts %ld best %ld (vertices in conflict %ld)\n", (long)(t - t0), it, conf, best, nl); tl = t; }
            if (difftime(t, t0) > maxsec) break;
        }
    }
    fprintf(stderr, "end it=%ld conflicts %ld best %ld\n", it, conf, best);
    char fn[256]; snprintf(fn, sizeof fn, "best_p%d_r%d_K%d_%s.txt", p, r, K, argv[5]);
    FILE *f = fopen(fn, "w"); for (long v = 0; v < n; v++) fprintf(f, "%d\n", bestcol[v]); fclose(f);
    printf("%s %ld\n", best == 0 ? "OK" : "FAIL", best);
    return 0;
}
