/* tabucol: K-colouring by tabu search (Hertz-de Werra, Galinier-Hao incremental table).
   stdin : n m K maxiter seed, then m edges "a b", then n initial colours (-1 = choose greedily)
   stdout: "OK <iters>" or "FAIL <best conflicts>", then the n colours of the best colouring.  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs;
static inline unsigned rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }
int main(void){
  int n, m, K; long long maxit; unsigned long long seed;
  if (scanf("%d %d %d %lld %llu", &n, &m, &K, &maxit, &seed) != 5) return 1;
  rs = seed * 2654435761ULL + 88172645463325252ULL;
  int *ea = malloc(sizeof(int) * m), *eb = malloc(sizeof(int) * m), *deg = calloc(n + 1, sizeof(int));
  for (int i = 0; i < m; i++){ if (scanf("%d %d", &ea[i], &eb[i]) != 2) return 1; deg[ea[i]]++; deg[eb[i]]++; }
  int *off = malloc(sizeof(int) * (n + 1)); off[0] = 0; for (int v = 0; v < n; v++) off[v + 1] = off[v] + deg[v];
  int *adj = malloc(sizeof(int) * (2 * (size_t)m + 1)), *pos = calloc(n, sizeof(int));
  for (int i = 0; i < m; i++){ adj[off[ea[i]] + pos[ea[i]]++] = eb[i]; adj[off[eb[i]] + pos[eb[i]]++] = ea[i]; }
  int *col = malloc(sizeof(int) * n);
  for (int v = 0; v < n; v++) if (scanf("%d", &col[v]) != 1) return 1;
  int *gam = calloc((size_t)n * K, sizeof(int));
  /* greedy for uncoloured vertices, in input order */
  for (int v = 0; v < n; v++) if (col[v] < 0) {
    int cnt[64] = {0};
    for (int j = off[v]; j < off[v + 1]; j++) if (col[adj[j]] >= 0) cnt[col[adj[j]]]++;
    int bc = 0; for (int c = 1; c < K; c++) if (cnt[c] < cnt[bc] || (cnt[c] == cnt[bc] && (rnd() & 1))) bc = c;
    col[v] = bc;
  }
  for (int v = 0; v < n; v++) for (int j = off[v]; j < off[v + 1]; j++) gam[(size_t)v * K + col[adj[j]]]++;
  long long f = 0; for (int v = 0; v < n; v++) f += gam[(size_t)v * K + col[v]]; f /= 2;
  int *inconf = malloc(sizeof(int) * n), *where = malloc(sizeof(int) * n), nc = 0;
  for (int v = 0; v < n; v++){ where[v] = -1; if (gam[(size_t)v * K + col[v]] > 0){ where[v] = nc; inconf[nc++] = v; } }
  long long *tabu = calloc((size_t)n * K, sizeof(long long));
  int *best = malloc(sizeof(int) * n); memcpy(best, col, sizeof(int) * n); long long bestf = f;
  long long it = 0;
  while (f > 0 && it < maxit){
    it++;
    int bv = -1, bcol = -1, bd = 1 << 30, ties = 0;
    for (int i = 0; i < nc; i++){
      int v = inconf[i]; int cv = col[v]; int gv = gam[(size_t)v * K + cv];
      for (int c = 0; c < K; c++){ if (c == cv) continue;
        int d = gam[(size_t)v * K + c] - gv;
        int istabu = tabu[(size_t)v * K + c] > it;
        if (istabu && !(f + d < bestf)) continue;
        if (d < bd){ bd = d; bv = v; bcol = c; ties = 1; }
        else if (d == bd){ ties++; if (rnd() % ties == 0){ bv = v; bcol = c; } }
      }
    }
    if (bv < 0){ /* everything tabu: random move among conflicting vertices */
      bv = inconf[rnd() % nc]; bcol = rnd() % K; if (bcol == col[bv]) bcol = (bcol + 1) % K;
      bd = gam[(size_t)bv * K + bcol] - gam[(size_t)bv * K + col[bv]];
    }
    int old = col[bv]; col[bv] = bcol; f += bd;
    for (int j = off[bv]; j < off[bv + 1]; j++){
      int w = adj[j];
      gam[(size_t)w * K + old]--; gam[(size_t)w * K + bcol]++;
      int inw = gam[(size_t)w * K + col[w]] > 0;
      if (inw && where[w] < 0){ where[w] = nc; inconf[nc++] = w; }
      else if (!inw && where[w] >= 0){ int p = where[w]; inconf[p] = inconf[--nc]; where[inconf[p]] = p; where[w] = -1; }
    }
    { int inv = gam[(size_t)bv * K + col[bv]] > 0;
      if (inv && where[bv] < 0){ where[bv] = nc; inconf[nc++] = bv; }
      else if (!inv && where[bv] >= 0){ int p = where[bv]; inconf[p] = inconf[--nc]; where[inconf[p]] = p; where[bv] = -1; } }
    tabu[(size_t)bv * K + old] = it + (rnd() % 10) + (long long)(0.6 * nc);
    if (f < bestf){ bestf = f; memcpy(best, col, sizeof(int) * n); }
  }
  if (f == 0) printf("OK %lld\n", it); else printf("FAIL %lld\n", bestf);
  for (int v = 0; v < n; v++) printf("%d\n", f == 0 ? col[v] : best[v]);
  return 0;
}
