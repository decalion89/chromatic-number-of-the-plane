/* tabucol.c -- tabu search for a proper k-colouring (Hertz & de Werra), starting from a given colouring.
   input file: n m k maxiter seed nfixed
               then nfixed lines "v c" (fixed vertices), then n initial colours (-1 = none), then m edges "a b"
   output (stdout): "OK" + n colours, or "FAIL best_conflicts"
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long rs = 88172645463325252ULL;
static inline unsigned long long xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r");
    if (!f) return 2;
    int n, m, k; long long maxit; unsigned long long seed; int nfix;
    if (fscanf(f, "%d %d %d %lld %llu %d", &n, &m, &k, &maxit, &seed, &nfix) != 6) return 3;
    rs ^= seed * 2654435761ULL + 1;
    char *fixed = calloc(n, 1);
    int *col = malloc(sizeof(int) * n);
    int *fv = malloc(sizeof(int) * (nfix + 1)), *fc = malloc(sizeof(int) * (nfix + 1));
    for (int i = 0; i < nfix; i++) { if (fscanf(f, "%d %d", &fv[i], &fc[i]) != 2) return 4; }
    for (int i = 0; i < n; i++) { if (fscanf(f, "%d", &col[i]) != 1) return 5; }
    int *ea = malloc(sizeof(int) * m), *eb = malloc(sizeof(int) * m);
    int *deg = calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) { if (fscanf(f, "%d %d", &ea[i], &eb[i]) != 2) return 6; deg[ea[i]]++; deg[eb[i]]++; }
    fclose(f);
    int *off = malloc(sizeof(int) * (n + 1)); off[0] = 0;
    for (int i = 0; i < n; i++) off[i + 1] = off[i] + deg[i];
    int *adj = malloc(sizeof(int) * 2 * m), *pos = calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) { adj[off[ea[i]] + pos[ea[i]]++] = eb[i]; adj[off[eb[i]] + pos[eb[i]]++] = ea[i]; }
    for (int i = 0; i < nfix; i++) { fixed[fv[i]] = 1; col[fv[i]] = fc[i]; }
    /* uncoloured vertices: greedy least-conflict colour */
    for (int v = 0; v < n; v++) if (col[v] < 0 || col[v] >= k) {
        int cnt[16] = {0};
        for (int e = off[v]; e < off[v + 1]; e++) { int w = adj[e]; if (col[w] >= 0 && col[w] < k) cnt[col[w]]++; }
        int best = 0; for (int c = 1; c < k; c++) if (cnt[c] < cnt[best] || (cnt[c] == cnt[best] && (xr() & 1))) best = c;
        col[v] = best;
    }
    int *g = calloc((size_t)n * k, sizeof(int));      /* g[v*k+c] = #neighbours of v with colour c */
    for (int v = 0; v < n; v++) for (int e = off[v]; e < off[v + 1]; e++) g[(size_t)v * k + col[adj[e]]]++;
    long long *tabu = calloc((size_t)n * k, sizeof(long long));
    int *cv = malloc(sizeof(int) * n), *cpos = malloc(sizeof(int) * n), ncv = 0;   /* conflicting vertices */
    for (int v = 0; v < n; v++) cpos[v] = -1;
    long long fcur = 0;
    for (int v = 0; v < n; v++) if (g[(size_t)v * k + col[v]] > 0) { fcur += g[(size_t)v * k + col[v]]; cpos[v] = ncv; cv[ncv++] = v; }
    fcur /= 2;
    long long fbest = fcur; int *bestcol = malloc(sizeof(int) * n); memcpy(bestcol, col, sizeof(int) * n);
    for (long long it = 0; it < maxit && fcur > 0; it++) {
        int bv = -1, bc = -1, bd = 1 << 30, nb = 0;
        for (int i = 0; i < ncv; i++) {
            int v = cv[i]; if (fixed[v]) continue;
            int cc = col[v]; int base = g[(size_t)v * k + cc];
            for (int c = 0; c < k; c++) if (c != cc) {
                int d = g[(size_t)v * k + c] - base;
                int ok = tabu[(size_t)v * k + c] <= it || fcur + d < fbest;
                if (!ok) continue;
                if (d < bd) { bd = d; bv = v; bc = c; nb = 1; }
                else if (d == bd) { nb++; if (xr() % nb == 0) { bv = v; bc = c; } }
            }
        }
        if (bv < 0) { /* all moves tabu: random move of a conflicting non-fixed vertex */
            int tries = 0;
            do { bv = cv[xr() % ncv]; tries++; } while (fixed[bv] && tries < 100);
            if (fixed[bv]) break;
            bc = (col[bv] + 1 + (int)(xr() % (k - 1))) % k;
            bd = g[(size_t)bv * k + bc] - g[(size_t)bv * k + col[bv]];
        }
        int old = col[bv];
        col[bv] = bc; fcur += bd;
        tabu[(size_t)bv * k + old] = it + (long long)(xr() % 10) + (long long)(0.6 * ncv);
        for (int e = off[bv]; e < off[bv + 1]; e++) {
            int w = adj[e];
            g[(size_t)w * k + old]--; g[(size_t)w * k + bc]++;
            int conf = g[(size_t)w * k + col[w]] > 0;
            if (conf && cpos[w] < 0) { cpos[w] = ncv; cv[ncv++] = w; }
            else if (!conf && cpos[w] >= 0) { int j = cpos[w]; cv[j] = cv[--ncv]; cpos[cv[j]] = j; cpos[w] = -1; }
        }
        {
            int conf = g[(size_t)bv * k + bc] > 0;
            if (conf && cpos[bv] < 0) { cpos[bv] = ncv; cv[ncv++] = bv; }
            else if (!conf && cpos[bv] >= 0) { int j = cpos[bv]; cv[j] = cv[--ncv]; cpos[cv[j]] = j; cpos[bv] = -1; }
        }
        if (fcur < fbest) { fbest = fcur; memcpy(bestcol, col, sizeof(int) * n); }
    }
    if (fbest == 0) {
        printf("OK\n");
        for (int v = 0; v < n; v++) printf("%d\n", bestcol[v]);
    } else printf("FAIL %lld\n", fbest);
    return 0;
}
