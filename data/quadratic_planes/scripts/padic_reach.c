/* padic_reach.c: the primes p = 3 (mod 4) that none of the certified fields Q(sqrt d) reaches.
 *
 * chi(Q_p^2) >= 4 as soon as one d below (each Q(sqrt d) carries a certified unit-distance graph with no
 * 3-colouring, data/quadratic_planes/q{d}.json) is a nonzero square mod p, since then Q(sqrt d) embeds in Q_p.
 * This program lists, in increasing order, the primes p = 3 (mod 4) up to LIMIT for which no d is a nonzero
 * square mod p (Jacobi symbols over a segmented sieve of Eratosthenes). padic_reach.py checks its list of d
 * against the data and its first answers again.
 *
 * build: cc -O2 -o padic_reach padic_reach.c      run: ./padic_reach LIMIT [COUNT]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static const int D[] = {11, 23, 35, 47, 59, 71, 95, 119, 131, 155, 179, 191, 239, 251, 263, 359, 431, 443, 455,
                        491, 599, 611, 791, 851, 911, 935, 959};
#define ND ((int)(sizeof D / sizeof D[0]))

static int jacobi(uint64_t a, uint64_t n) {   /* n odd, n > 0 */
  int t = 1;
  a %= n;
  while (a) {
    while ((a & 1) == 0) { a >>= 1; uint64_t r = n & 7; if (r == 3 || r == 5) t = -t; }
    uint64_t tmp = a; a = n; n = tmp;
    if ((a & 3) == 3 && (n & 3) == 3) t = -t;
    a %= n;
  }
  return n == 1 ? t : 0;
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s LIMIT [COUNT]\n", argv[0]); return 2; }
  uint64_t limit = strtoull(argv[1], 0, 10);
  int want = argc > 2 ? atoi(argv[2]) : 5, found = 0;
  uint64_t seglen = (uint64_t)1 << 24, root = 1, count = 0;
  while (root * root <= limit) root++;
  char *composite = calloc(root + 1, 1);
  uint64_t *primes = malloc(sizeof(uint64_t) * (root + 1));
  int np = 0;
  for (uint64_t i = 2; i <= root; i++)
    if (!composite[i]) { primes[np++] = i; for (uint64_t j = i * i; j <= root; j += i) composite[j] = 1; }
  char *seg = malloc(seglen);
  printf("d = %d fields; primes p = 3 (mod 4) up to %llu with no d a nonzero square mod p:\n", ND,
         (unsigned long long)limit);
  for (uint64_t lo = 2; lo <= limit && found < want; lo += seglen) {
    uint64_t hi = lo + seglen - 1;
    if (hi > limit) hi = limit;
    memset(seg, 0, hi - lo + 1);
    for (int k = 0; k < np && primes[k] * primes[k] <= hi; k++) {
      uint64_t p = primes[k], start = (lo + p - 1) / p * p;
      if (start < p * p) start = p * p;
      for (uint64_t j = start; j <= hi; j += p) seg[j - lo] = 1;
    }
    for (uint64_t x = lo; x <= hi; x++) {
      if (seg[x - lo] || (x & 3) != 3) continue;
      count++;
      int reached = 0;
      for (int i = 0; i < ND && !reached; i++) reached = jacobi((uint64_t)D[i], x) == 1;
      if (!reached) {
        printf("p = %llu (number %llu among the primes = 3 mod 4)\n", (unsigned long long)x, (unsigned long long)count);
        fflush(stdout);
        if (++found >= want) break;
      }
    }
  }
  printf("searched to %llu\n", (unsigned long long)limit);
  return 0;
}
