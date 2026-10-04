/* Exhaustive: count maps f: V -> Z/8 with f(0) = 0 that are homomorphisms to K_{8/3} (|f(a)-f(b)| mod 8 in {3,4,5}),
   for M, H7 and the q31 graph (edge lists in the labels of the witness files), and for every H7 - v. */
#include <stdio.h>
static int ok(int a, int b) { int t = ((b - a) % 8 + 8) % 8; return t >= 3 && t <= 5; }
static long count(int ne, int E[][2], int skip) {
    long cnt = 0; int f[9];
    for (long m = 0; m < (1L << 24); m++) {          /* 8^8 maps of vertices 1..8, f(0) = 0 */
        f[0] = 0; long t = m;
        for (int v = 1; v < 9; v++) { f[v] = t & 7; t >>= 3; }
        int good = 1;
        for (int e = 0; e < ne && good; e++) {
            if (E[e][0] == skip || E[e][1] == skip) continue;
            if (!ok(f[E[e][0]], f[E[e][1]])) good = 0;
        }
        cnt += good;
    }
    return cnt;
}
int main(void) {
    int H7[13][2] = {{0,1},{0,4},{0,7},{1,2},{1,5},{2,3},{2,6},{3,4},{3,8},{4,5},{5,6},{6,7},{7,8}};
    int M7[12][2] = {{0,1},{0,4},{0,7},{1,2},{2,3},{2,6},{3,4},{3,8},{4,5},{5,6},{6,7},{7,8}};   /* H7 - P1P5 */
    int Q31[14][2] = {{0,1},{0,3},{0,6},{1,4},{1,8},{2,4},{2,6},{2,7},{3,4},{3,7},{4,5},{5,6},{5,8},{7,8}};
    printf("homomorphisms to K_8/3 with f(0)=0: M %ld, H7 %ld, q31 graph %ld\n", count(12, M7, -1), count(13, H7, -1), count(14, Q31, -1));
    printf("H7 - v (v = 1..8; v = 0 would need another fixed vertex):");
    for (int v = 1; v < 9; v++) printf(" %ld", count(13, H7, v));
    printf("\n");
    return 0;
}
