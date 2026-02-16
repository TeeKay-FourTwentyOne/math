/*
 * Verify min P(b)/Q(b) > 1 for all achievable B-profiles at p=31.
 *
 * Enumerates all C(31,15) = 300,540,195 D12 subsets of Z_31.
 * For each D12, computes B-profile (B(1),...,B(15)) where B(d) = |D12 ∩ (D12+d)|.
 * All profiles lie on hyperplane Σ B(d) = S = 15*14/2 = 105.
 *
 * Stores joint distribution in hash table, marginals in arrays.
 * Then computes P(b)/Q(b) for each profile and reports minimum.
 *
 * If min P/Q > 1: validates that hyperplane conditioning bound is a valid
 * lower bound on E[N], proving existence for primes up to p~227.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define P 31
#define K 15          /* (P-1)/2 = |D12| */
#define R 15          /* number of B-components */
#define S 105         /* K*(K-1)/2 = hyperplane sum */
#define MASK_P ((1u << P) - 1)

/* Hash table for B-profile joint distribution */
#define HT_BITS 24
#define HT_SIZE (1u << HT_BITS)
#define HT_MASK (HT_SIZE - 1)
#define EMPTY_KEY 0xFFFFFFFFFFFFFFFFULL

typedef struct {
    uint64_t key;
    uint32_t count;
} ht_entry_t;

static ht_entry_t *ht;
static uint32_t ht_used = 0;

/* Marginal counts: marginal[j][b] = # of D12 with B(j+1) = b */
static uint64_t marginal[R][K];  /* B(d) ∈ {0,...,K-1} */

static inline uint32_t ht_hash(uint64_t key) {
    /* murmur3 finalizer */
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return (uint32_t)(key & HT_MASK);
}

static inline void ht_insert(uint64_t key) {
    uint32_t idx = ht_hash(key);
    while (1) {
        if (ht[idx].key == EMPTY_KEY) {
            ht[idx].key = key;
            ht[idx].count = 1;
            ht_used++;
            return;
        }
        if (ht[idx].key == key) {
            ht[idx].count++;
            return;
        }
        idx = (idx + 1) & HT_MASK;
    }
}

static inline uint32_t rot(uint32_t m, int d) {
    return ((m << d) | (m >> (P - d))) & MASK_P;
}

static inline uint32_t next_comb(uint32_t v) {
    uint32_t t = v | (v - 1);
    return (t + 1) | (((~t & -(~t)) - 1) >> (__builtin_ctz(v) + 1));
}

int main(void) {
    printf("=== Min P/Q Verification at p=%d ===\n", P);
    printf("C(%d,%d) = 300,540,195 subsets to enumerate\n", P, K);
    printf("R=%d components, hyperplane sum S=%d\n", R, S);
    printf("Hash table: %u slots (%.0f MB)\n\n",
           HT_SIZE, (double)HT_SIZE * sizeof(ht_entry_t) / 1e6);

    /* Allocate hash table */
    ht = (ht_entry_t *)malloc(HT_SIZE * sizeof(ht_entry_t));
    if (!ht) { fprintf(stderr, "OOM for hash table\n"); return 1; }
    for (uint32_t i = 0; i < HT_SIZE; i++) {
        ht[i].key = EMPTY_KEY;
        ht[i].count = 0;
    }
    memset(marginal, 0, sizeof(marginal));

    /* Precompute nothing — all done inline */
    time_t t0 = time(NULL);
    uint64_t total = 0;
    uint32_t d12 = (1u << K) - 1;  /* first K-subset: bits 0..14 */

    while (1) {
        /* Compute B-profile and pack into 60-bit key */
        uint64_t key = 0;
        for (int j = 0; j < R; j++) {
            int d = j + 1;
            uint32_t rd = rot(d12, d);
            int b = __builtin_popcount(d12 & rd);
            key |= ((uint64_t)b) << (j * 4);
            marginal[j][b]++;
        }

        ht_insert(key);
        total++;

        if ((total & 0xFFFFFF) == 0) {
            time_t now = time(NULL);
            long elapsed = now - t0;
            double pct = total * 100.0 / 300540195.0;
            double rate = elapsed > 0 ? total / (double)elapsed : 0;
            long eta = rate > 0 ? (long)((300540195 - total) / rate) : 0;
            printf("  %9llu / 300540195 (%5.1f%%) [%lds, ETA %lds] ht=%u\n",
                   (unsigned long long)total, pct, elapsed, eta, ht_used);
            fflush(stdout);
        }

        d12 = next_comb(d12);
        if (d12 > MASK_P) break;
    }

    time_t t1 = time(NULL);
    printf("\nEnumeration complete: %llu subsets in %lds\n",
           (unsigned long long)total, t1 - t0);
    printf("Distinct B-profiles: %u\n", ht_used);
    printf("Hash table occupancy: %.1f%%\n\n", 100.0 * ht_used / HT_SIZE);

    if (ht_used > HT_SIZE * 3 / 4) {
        fprintf(stderr, "WARNING: hash table >75%% full!\n");
    }

    /* Verify hyperplane constraint */
    int hp_ok = 1;
    for (uint32_t i = 0; i < HT_SIZE; i++) {
        if (ht[i].key == EMPTY_KEY) continue;
        int sum = 0;
        for (int j = 0; j < R; j++)
            sum += (ht[i].key >> (j * 4)) & 0xF;
        if (sum != S) {
            printf("ERROR: profile sum=%d != %d\n", sum, S);
            hp_ok = 0;
            break;
        }
    }
    printf("Hyperplane check (all sums = %d): %s\n", S, hp_ok ? "PASS" : "FAIL");

    /* Verify marginals are identical (by cyclic symmetry of Z_p) */
    int marg_identical = 1;
    for (int j = 1; j < R; j++) {
        for (int b = 0; b < K; b++) {
            if (marginal[j][b] != marginal[0][b]) {
                marg_identical = 0;
                printf("Marginal mismatch: j=%d b=%d: %llu vs %llu\n",
                       j, b, (unsigned long long)marginal[j][b],
                       (unsigned long long)marginal[0][b]);
                break;
            }
        }
        if (!marg_identical) break;
    }
    printf("Marginals identical: %s\n", marg_identical ? "YES" : "NO");

    /* Print marginal distribution */
    printf("\nMarginal PMF (same for all d):\n");
    for (int b = 0; b < K; b++) {
        if (marginal[0][b] > 0)
            printf("  B(d)=%2d: %12llu (%.6e)\n",
                   b, (unsigned long long)marginal[0][b],
                   (double)marginal[0][b] / total);
    }

    /* Compute min P/Q ratio */
    printf("\n=== P/Q RATIO ANALYSIS ===\n");
    double min_ratio = 1e30, max_ratio = 0;
    uint64_t below_1 = 0;
    double sum_log_ratio = 0;
    uint64_t n_profiles = 0;
    uint64_t min_key = 0;
    uint32_t min_count = 0;

    for (uint32_t i = 0; i < HT_SIZE; i++) {
        if (ht[i].key == EMPTY_KEY) continue;

        /* P(b) = count / total */
        double log_P = log((double)ht[i].count) - log((double)total);

        /* Q(b) = Π_j marginal[j][b_j] / total */
        double log_Q = 0;
        for (int j = 0; j < R; j++) {
            int b = (ht[i].key >> (j * 4)) & 0xF;
            log_Q += log((double)marginal[j][b]) - log((double)total);
        }

        double log_ratio = log_P - log_Q;
        double ratio = exp(log_ratio);

        if (ratio < min_ratio) {
            min_ratio = ratio;
            min_key = ht[i].key;
            min_count = ht[i].count;
        }
        if (ratio > max_ratio) max_ratio = ratio;
        if (ratio < 1.0) below_1++;
        sum_log_ratio += log_ratio;
        n_profiles++;
    }

    printf("min P/Q = %.10f\n", min_ratio);
    printf("max P/Q = %.6e\n", max_ratio);
    printf("geometric mean P/Q = %.6f\n", exp(sum_log_ratio / n_profiles));
    printf("Profiles with P/Q < 1: %llu / %llu\n",
           (unsigned long long)below_1, (unsigned long long)n_profiles);
    printf("\n*** min P/Q > 1: %s ***\n", min_ratio > 1.0 ? "YES" : "NO");

    if (min_ratio > 1.0) {
        printf("log2(min P/Q) = %.4f\n", log(min_ratio) / log(2.0));
        printf("\nThis validates the hyperplane conditioning bound!\n");
        printf("E[N] >= 2^margin for all D11 with positive margin.\n");
    }

    /* Print min-ratio profile */
    printf("\nMin-ratio profile (count=%u, P=%.3e):\n  B = [",
           min_count, (double)min_count / total);
    for (int j = 0; j < R; j++) {
        int b = (min_key >> (j * 4)) & 0xF;
        printf("%d%s", b, j < R - 1 ? ", " : "");
    }
    printf("]\n");

    /* Print histogram of log2(P/Q) */
    printf("\nlog2(P/Q) histogram:\n");
    int hist[100];
    memset(hist, 0, sizeof(hist));
    for (uint32_t i = 0; i < HT_SIZE; i++) {
        if (ht[i].key == EMPTY_KEY) continue;
        double log_P = log((double)ht[i].count) - log((double)total);
        double log_Q = 0;
        for (int j = 0; j < R; j++) {
            int b = (ht[i].key >> (j * 4)) & 0xF;
            log_Q += log((double)marginal[j][b]) - log((double)total);
        }
        double log2_ratio = (log_P - log_Q) / log(2.0);
        int bin = (int)(log2_ratio + 0.5);
        if (bin >= 0 && bin < 100) hist[bin]++;
    }
    for (int b = 0; b < 100; b++) {
        if (hist[b] > 0) printf("  [%2d,%2d): %d\n", b, b + 1, hist[b]);
    }

    /* Write JSON results */
    FILE *fout = fopen("ramsey-book-graphs/min_pq_p31_results.json", "w");
    if (fout) {
        fprintf(fout, "{\n");
        fprintf(fout, "  \"p\": %d,\n", P);
        fprintf(fout, "  \"k\": %d,\n", K);
        fprintf(fout, "  \"total_subsets\": %llu,\n", (unsigned long long)total);
        fprintf(fout, "  \"distinct_profiles\": %u,\n", ht_used);
        fprintf(fout, "  \"min_PQ\": %.10f,\n", min_ratio);
        fprintf(fout, "  \"max_PQ\": %.6e,\n", max_ratio);
        fprintf(fout, "  \"below_1\": %llu,\n", (unsigned long long)below_1);
        fprintf(fout, "  \"min_PQ_gt_1\": %s,\n", min_ratio > 1.0 ? "true" : "false");
        fprintf(fout, "  \"log2_min_PQ\": %.6f\n", log(min_ratio) / log(2.0));
        fprintf(fout, "}\n");
        fclose(fout);
        printf("\nResults written to ramsey-book-graphs/min_pq_p31_results.json\n");
    }

    free(ht);
    return 0;
}
