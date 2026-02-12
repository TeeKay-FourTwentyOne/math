/*
 * Exact N(D11) enumeration at p=31 for second moment ratio.
 *
 * For each symmetric D11 orbit, count valid D12 subsets.
 * D12 subset of {0,...,30}, |D12| = 15 = (p-1)/2.
 * Total D12 candidates: C(31,15) = 300,540,195.
 *
 * Constraints (VALIDATED against ramsey_core.verify_construction):
 *   V1V1 red  (d in D11): A(d) + B(d) <= 14        [n-2]
 *   V1V1 blue (d in D22): A(d) + B(d) <= 17        [n+1, usually slack]
 *   V2V2 red  (d in D22): C(d) + B(p-d) <= 14      [n-2]
 *   V2V2 blue (d in D11): C(d) + B(p-d) <= 13      [n-3, TIGHTEST]
 *   V1V2 red  (d in D12): X(d) <= 14               [n-2]
 *   V1V2 blue (d not in D12): X(d) <= 15            [n-1]
 *
 * where B(d) = Delta(D12,D12,d), X(d) = Sigma(D11,D12,d) + Delta(D12,D22,d)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define P 31
#define N_BOOK 16
#define HALF_P 15
#define K_PAIRS 8
#define D12_SIZE 15
#define MASK_P ((1u << P) - 1)
#define PRIM_ROOT 3

static inline uint32_t rot(uint32_t m, int d) {
    d %= P;
    if (d == 0) return m;
    return ((m << d) | (m >> (P - d))) & MASK_P;
}

static inline uint32_t rev(uint32_t m) {
    uint32_t r = m & 1;
    for (int i = 1; i < P; i++)
        if (m & (1u << i)) r |= 1u << (P - i);
    return r;
}

static inline uint32_t next_comb(uint32_t v) {
    uint32_t t = v | (v - 1);
    return (t + 1) | (((~t & -(~t)) - 1) >> (__builtin_ctz(v) + 1));
}

static uint32_t mul_set(uint32_t m, int g) {
    uint32_t r = 0;
    for (int i = 0; i < P; i++)
        if (m & (1u << i)) r |= 1u << ((i * g) % P);
    return r;
}

int generate_orbits(uint32_t *reps, int *sizes) {
    int p0[HALF_P], p1[HALF_P];
    for (int i = 0; i < HALF_P; i++) {
        p0[i] = i + 1;
        p1[i] = P - (i + 1);
    }

    int n_d11 = 0, n_orb = 0;
    uint32_t *all_d11 = malloc(6500 * sizeof(uint32_t));
    uint32_t *canon = malloc(6500 * sizeof(uint32_t));

    uint32_t ps = (1u << K_PAIRS) - 1;
    uint32_t ps_limit = 1u << HALF_P;

    while (ps < ps_limit) {
        uint32_t d11 = 0;
        for (int i = 0; i < HALF_P; i++)
            if (ps & (1u << i)) {
                d11 |= 1u << p0[i];
                d11 |= 1u << p1[i];
            }
        all_d11[n_d11] = d11;

        uint32_t c = d11, cur = d11;
        for (int j = 0; j < HALF_P - 1; j++) {
            cur = mul_set(cur, PRIM_ROOT);
            if (cur < c) c = cur;
        }
        canon[n_d11] = c;
        n_d11++;

        ps = next_comb(ps);
        if (__builtin_popcount(ps) != K_PAIRS) break;
    }
    printf("Total D11: %d\n", n_d11);

    uint8_t *used = calloc(n_d11, 1);
    for (int i = 0; i < n_d11; i++) {
        if (used[i]) continue;
        reps[n_orb] = all_d11[i];
        uint32_t cur = all_d11[i];
        int sz = 1;
        for (int j = 0; j < HALF_P - 1; j++) {
            cur = mul_set(cur, PRIM_ROOT);
            if (cur == all_d11[i]) break;
            sz++;
        }
        sizes[n_orb] = sz;
        for (int j = i + 1; j < n_d11; j++)
            if (canon[j] == canon[i]) used[j] = 1;
        n_orb++;
    }

    free(all_d11); free(canon); free(used);
    printf("Orbits: %d\n", n_orb);
    return n_orb;
}

void validate_solution(uint32_t d11, uint32_t d12) {
    uint32_t d22 = MASK_P & ~d11 & ~1u;
    printf("Validate: |D11|=%d |D12|=%d |D22|=%d\n",
           __builtin_popcount(d11), __builtin_popcount(d12), __builtin_popcount(d22));
    int ok = 1;
    for (int d = 1; d < P; d++) {
        int B = __builtin_popcount(d12 & rot(d12, d));
        int Bpd = __builtin_popcount(d12 & rot(d12, P - d));
        if (d11 & (1u << d)) {
            int A = __builtin_popcount(d11 & rot(d11, d));
            int C = __builtin_popcount(d22 & rot(d22, d));
            if (A + B > N_BOOK - 2) { printf("  FAIL V1V1r d=%d: %d\n", d, A+B); ok=0; }
            if (C + Bpd > N_BOOK - 3) { printf("  FAIL V2V2b d=%d: %d\n", d, C+Bpd); ok=0; }
        } else {
            int A = __builtin_popcount(d11 & rot(d11, d));
            int C = __builtin_popcount(d22 & rot(d22, d));
            if (C + Bpd > N_BOOK - 2) { printf("  FAIL V2V2r d=%d: %d\n", d, C+Bpd); ok=0; }
            if (A + B > N_BOOK + 1) { printf("  FAIL V1V1b d=%d: %d\n", d, A+B); ok=0; }
        }
    }
    uint32_t d12r = rev(d12);
    uint32_t d22v = MASK_P & ~d11 & ~1u;
    for (int d = 0; d < P; d++) {
        int sig = __builtin_popcount(d11 & rot(d12r, d));
        int del = __builtin_popcount(d12 & rot(d22v, d));
        int X = sig + del;
        if (d12 & (1u << d)) {
            if (X > N_BOOK - 2) { printf("  FAIL V1V2r d=%d: %d\n", d, X); ok=0; }
        } else {
            if (X > N_BOOK - 1) { printf("  FAIL V1V2b d=%d: %d\n", d, X); ok=0; }
        }
    }
    printf("Result: %s\n\n", ok ? "PASS" : "FAIL");
}

int main(void) {
    /* Validate known solution */
    {
        int d11e[] = {2,3,4,9,10,12,13,14,17,18,19,21,22,27,28,29};
        int d12e[] = {0,2,4,5,7,9,10,11,12,13,20,24,27,28,30};
        uint32_t d11 = 0, d12 = 0;
        for (int i = 0; i < 16; i++) d11 |= 1u << d11e[i];
        for (int i = 0; i < 15; i++) d12 |= 1u << d12e[i];
        validate_solution(d11, d12);
    }

    uint32_t reps[500];
    int sizes[500];
    int n_orb = generate_orbits(reps, sizes);

    int total_d11 = 0;
    for (int i = 0; i < n_orb; i++) total_d11 += sizes[i];
    printf("Sum orbit sizes: %d (expect 6435)\n\n", total_d11);

    double sum_N = 0, sum_N2 = 0;
    long long total_weight = 0;
    FILE *fout = fopen("exact_N_p31_results.jsonl", "w");
    time_t t0 = time(NULL);

    for (int orb = 0; orb < n_orb; orb++) {
        uint32_t d11 = reps[orb];
        uint32_t d22 = MASK_P & ~d11 & ~1u;

        int A[P], Cv[P];
        int maxB_v1r[P], maxB_v2b[P];
        int maxB_v2r[P], maxB_v1b[P];
        int impossible = 0;

        for (int d = 1; d < P; d++) {
            A[d] = __builtin_popcount(d11 & rot(d11, d));
            Cv[d] = __builtin_popcount(d22 & rot(d22, d));
            if (d11 & (1u << d)) {
                maxB_v1r[d] = (N_BOOK - 2) - A[d];
                maxB_v2b[d] = (N_BOOK - 3) - Cv[d];
                if (maxB_v1r[d] < 0 || maxB_v2b[d] < 0) { impossible = 1; break; }
            } else {
                maxB_v2r[d] = (N_BOOK - 2) - Cv[d];
                maxB_v1b[d] = (N_BOOK + 1) - A[d];
                if (maxB_v2r[d] < 0) { impossible = 1; break; }
            }
        }

        if (impossible) {
            fprintf(fout, "{\"orbit\":%d,\"size\":%d,\"N\":0}\n", orb, sizes[orb]);
            fflush(fout);
            total_weight += sizes[orb];
            if ((orb+1) % 50 == 0)
                printf("Orbit %d/%d: N=0 (impossible) [%lds]\n",
                       orb+1, n_orb, time(NULL)-t0);
            continue;
        }

        /* Order d-values by tightness for early exit */
        int d_order[P-1], d_thr[P-1];
        for (int d = 1; d < P; d++) {
            d_order[d-1] = d;
            if (d11 & (1u << d))
                d_thr[d-1] = maxB_v1r[d] < maxB_v2b[d] ? maxB_v1r[d] : maxB_v2b[d];
            else
                d_thr[d-1] = maxB_v2r[d] < maxB_v1b[d] ? maxB_v2r[d] : maxB_v1b[d];
        }
        for (int i = 1; i < P-1; i++) {
            int td = d_order[i], tt = d_thr[i];
            int j = i;
            while (j > 0 && d_thr[j-1] > tt) {
                d_order[j] = d_order[j-1];
                d_thr[j] = d_thr[j-1];
                j--;
            }
            d_order[j] = td;
            d_thr[j] = tt;
        }

        uint32_t d22_rot[P];
        for (int d = 0; d < P; d++) d22_rot[d] = rot(d22, d);

        long long count = 0;
        uint32_t d12 = (1u << D12_SIZE) - 1;

        while (1) {
            int valid = 1;
            for (int i = 0; i < P-1 && valid; i++) {
                int d = d_order[i];
                int B_d = __builtin_popcount(d12 & rot(d12, d));
                int B_pd = __builtin_popcount(d12 & rot(d12, P - d));
                if (d11 & (1u << d)) {
                    if (B_d > maxB_v1r[d] || B_pd > maxB_v2b[d]) valid = 0;
                } else {
                    if (B_pd > maxB_v2r[d] || B_d > maxB_v1b[d]) valid = 0;
                }
            }

            if (valid) {
                uint32_t d12r = rev(d12);
                for (int d = 0; d < P && valid; d++) {
                    int sig = __builtin_popcount(d11 & rot(d12r, d));
                    int del = __builtin_popcount(d12 & d22_rot[d]);
                    int X = sig + del;
                    if (d12 & (1u << d)) {
                        if (X > N_BOOK - 2) valid = 0;
                    } else {
                        if (X > N_BOOK - 1) valid = 0;
                    }
                }
            }

            if (valid) count++;

            d12 = next_comb(d12);
            if (d12 > MASK_P) break;
        }

        fprintf(fout, "{\"orbit\":%d,\"size\":%d,\"N\":%lld}\n", orb, sizes[orb], count);
        fflush(fout);

        sum_N += (double)count * sizes[orb];
        sum_N2 += (double)count * count * sizes[orb];
        total_weight += sizes[orb];

        time_t now = time(NULL);
        long elapsed = now - t0;
        double frac = (double)(orb+1) / n_orb;
        long eta = frac > 0.01 ? (long)(elapsed / frac * (1 - frac)) : 0;
        printf("Orbit %d/%d: N=%lld (size=%d) [%lds, ETA %lds]\n",
               orb+1, n_orb, count, sizes[orb], elapsed, eta);
    }

    fclose(fout);

    double E_N = sum_N / total_weight;
    double E_N2 = sum_N2 / total_weight;
    double ratio = (E_N > 0) ? E_N2 / (E_N * E_N) : 0;

    printf("\n=== RESULTS ===\n");
    printf("Orbits: %d, Total D11: %lld\n", n_orb, total_weight);
    printf("E[N]   = %.6f\n", E_N);
    printf("E[N^2] = %.6f\n", E_N2);
    printf("E[N^2]/E[N]^2 = %.6f\n", ratio);

    return 0;
}
