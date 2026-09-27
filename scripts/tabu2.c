/* tabu2: proper K-colouring first (hard edges), then minimise soft conflicts WITHOUT ever
   breaking a hard edge (moves restricted to colours free of hard neighbours).
   stdin : n mh ms K it1 it2 seed ; mh hard edges "a b" ; ms soft edges "a b" ; n initial colours
   stdout: "OK <soft>" (proper; soft = alike soft edges of the best colouring) or "FAIL <hard>"; then colours */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs;
static inline unsigned rnd(void){ rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }
typedef struct { int *off, *adj; } csr;
static csr build(int n, int m, int *ea, int *eb){
  csr g; int *deg = calloc(n + 1, sizeof(int)); for (int i = 0; i < m; i++){ deg[ea[i]]++; deg[eb[i]]++; }
  g.off = malloc(sizeof(int) * (n + 1)); g.off[0] = 0; for (int v = 0; v < n; v++) g.off[v + 1] = g.off[v] + deg[v];
  g.adj = malloc(sizeof(int) * (2 * (size_t)m + 1)); int *pos = calloc(n, sizeof(int));
  for (int i = 0; i < m; i++){ g.adj[g.off[ea[i]] + pos[ea[i]]++] = eb[i]; g.adj[g.off[eb[i]] + pos[eb[i]]++] = ea[i]; }
  free(deg); free(pos); return g;
}
int main(void){
  int n, mh, ms, K; long long it1, it2; unsigned long long seed;
  if (scanf("%d %d %d %d %lld %lld %llu", &n, &mh, &ms, &K, &it1, &it2, &seed) != 7) return 1;
  rs = seed * 2654435761ULL + 88172645463325252ULL;
  int *ha = malloc(sizeof(int) * (mh + 1)), *hb = malloc(sizeof(int) * (mh + 1)), *sa = malloc(sizeof(int) * (ms + 1)), *sb = malloc(sizeof(int) * (ms + 1));
  for (int i = 0; i < mh; i++) if (scanf("%d %d", &ha[i], &hb[i]) != 2) return 1;
  for (int i = 0; i < ms; i++) if (scanf("%d %d", &sa[i], &sb[i]) != 2) return 1;
  csr H = build(n, mh, ha, hb), S = build(n, ms, sa, sb);
  int *col = malloc(sizeof(int) * n); for (int v = 0; v < n; v++) if (scanf("%d", &col[v]) != 1) return 1;
  for (int v = 0; v < n; v++) if (col[v] < 0){
    int cnt[64] = {0}; for (int j = H.off[v]; j < H.off[v + 1]; j++) if (col[H.adj[j]] >= 0) cnt[col[H.adj[j]]]++;
    int bc = 0; for (int c = 1; c < K; c++) if (cnt[c] < cnt[bc] || (cnt[c] == cnt[bc] && (rnd() & 1))) bc = c; col[v] = bc; }
  int *gh = calloc((size_t)n * K, sizeof(int)), *gs = calloc((size_t)n * K, sizeof(int));
  for (int v = 0; v < n; v++){ for (int j = H.off[v]; j < H.off[v + 1]; j++) gh[(size_t)v * K + col[H.adj[j]]]++;
                               for (int j = S.off[v]; j < S.off[v + 1]; j++) gs[(size_t)v * K + col[S.adj[j]]]++; }
  long long fh = 0, fs = 0; for (int v = 0; v < n; v++){ fh += gh[(size_t)v * K + col[v]]; fs += gs[(size_t)v * K + col[v]]; } fh /= 2; fs /= 2;
  long long *tabu = calloc((size_t)n * K, sizeof(long long));
  int *inc = malloc(sizeof(int) * n), *wh = malloc(sizeof(int) * n), nc = 0;
  #define MOVE(v, c) do { int _o = col[v]; col[v] = (c); \
      for (int j = H.off[v]; j < H.off[v + 1]; j++){ int w = H.adj[j]; gh[(size_t)w * K + _o]--; gh[(size_t)w * K + (c)]++; } \
      for (int j = S.off[v]; j < S.off[v + 1]; j++){ int w = S.adj[j]; gs[(size_t)w * K + _o]--; gs[(size_t)w * K + (c)]++; } } while (0)
  /* phase 1: hard conflicts only */
  long long it = 0;
  for (int v = 0; v < n; v++){ wh[v] = -1; if (gh[(size_t)v * K + col[v]] > 0){ wh[v] = nc; inc[nc++] = v; } }
  while (fh > 0 && it < it1){
    it++; int bv = -1, bc = -1, bd = 1 << 30, ties = 0;
    for (int i = 0; i < nc; i++){ int v = inc[i], cv = col[v], g0 = gh[(size_t)v * K + cv];
      for (int c = 0; c < K; c++){ if (c == cv) continue; int d = gh[(size_t)v * K + c] - g0;
        if (tabu[(size_t)v * K + c] > it && !(fh + d < 1)) continue;
        if (d < bd){ bd = d; bv = v; bc = c; ties = 1; } else if (d == bd){ ties++; if (rnd() % ties == 0){ bv = v; bc = c; } } } }
    if (bv < 0){ bv = inc[rnd() % nc]; bc = rnd() % K; if (bc == col[bv]) bc = (bc + 1) % K; bd = gh[(size_t)bv * K + bc] - gh[(size_t)bv * K + col[bv]]; }
    int old = col[bv];
    fs += gs[(size_t)bv * K + bc] - gs[(size_t)bv * K + old];
    MOVE(bv, bc); fh += bd;
    for (int j = H.off[bv]; j <= H.off[bv + 1]; j++){ int w = (j < H.off[bv + 1]) ? H.adj[j] : bv; int in = gh[(size_t)w * K + col[w]] > 0;
      if (in && wh[w] < 0){ wh[w] = nc; inc[nc++] = w; } else if (!in && wh[w] >= 0){ int p = wh[w]; inc[p] = inc[--nc]; wh[inc[p]] = p; wh[w] = -1; } }
    tabu[(size_t)bv * K + old] = it + (rnd() % 10) + (long long)(0.6 * nc);
  }
  if (fh > 0){ printf("FAIL %lld\n", fh); for (int v = 0; v < n; v++) printf("%d\n", col[v]); return 0; }
  /* phase 2: minimise W*hard + soft over hard-conflicting and soft-edge vertices; keep the best PROPER colouring */
  memset(tabu, 0, sizeof(long long) * (size_t)n * K);
  int *best = malloc(sizeof(int) * n); memcpy(best, col, sizeof(int) * n); long long bestfs = fs;
  int *sv = malloc(sizeof(int) * n), nsv = 0; for (int v = 0; v < n; v++) if (S.off[v + 1] > S.off[v]) sv[nsv++] = v;
  const int W = 8;
  for (long long t = 1; t <= it2 && bestfs > 0; t++){
    int bv = -1, bc = -1, bd = 1 << 30, ties = 0;
    for (int pass = 0; pass < 2; pass++){
      int cnt = pass == 0 ? nsv : nc;
      for (int ii = 0; ii < cnt; ii++){ int v = pass == 0 ? sv[ii] : inc[ii]; int cv = col[v];
        int h0 = gh[(size_t)v * K + cv], s0 = gs[(size_t)v * K + cv]; if (h0 == 0 && s0 == 0) continue;
        for (int c = 0; c < K; c++){ if (c == cv) continue;
          int dh = gh[(size_t)v * K + c] - h0, ds = gs[(size_t)v * K + c] - s0, d = W * dh + ds;
          int asp = (fh + dh == 0) && (fs + ds < bestfs);
          if (tabu[(size_t)v * K + c] > t && !asp) continue;
          if (d < bd){ bd = d; bv = v; bc = c; ties = 1; } else if (d == bd){ ties++; if (rnd() % ties == 0){ bv = v; bc = c; } } } } }
    if (bv < 0) break;
    int old = col[bv], dh = gh[(size_t)bv * K + bc] - gh[(size_t)bv * K + old], ds = gs[(size_t)bv * K + bc] - gs[(size_t)bv * K + old];
    MOVE(bv, bc); fh += dh; fs += ds;
    for (int j = H.off[bv]; j <= H.off[bv + 1]; j++){ int w = (j < H.off[bv + 1]) ? H.adj[j] : bv; int in = gh[(size_t)w * K + col[w]] > 0;
      if (in && wh[w] < 0){ wh[w] = nc; inc[nc++] = w; } else if (!in && wh[w] >= 0){ int p = wh[w]; inc[p] = inc[--nc]; wh[inc[p]] = p; wh[w] = -1; } }
    tabu[(size_t)bv * K + old] = t + (rnd() % 10) + 7 + nc;
    if (fh == 0 && fs < bestfs){ bestfs = fs; memcpy(best, col, sizeof(int) * n); }
  }
  printf("OK %lld\n", bestfs); for (int v = 0; v < n; v++) printf("%d\n", best[v]);
  return 0;
}
