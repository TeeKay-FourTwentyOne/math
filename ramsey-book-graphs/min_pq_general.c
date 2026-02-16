/*
 * Verify min P(b)/Q(b) > 1 for B-profiles of uniform k-subsets of Z_p.
 * GENERAL VERSION: works for any prime p (passed as argument).
 *
 * For p <= 31: B-profile packed into single uint64_t (4 bits × 15 = 60 bits)
 * For p <= 43: B-profile packed into two uint64_t (5 bits × 21 = 105 bits)
 *
 * Usage: ./min_pq_general <p>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Parameters set at runtime */
static int P_VAL, K_VAL, R_VAL;
static int S_VAL;   /* K*(K-1)/2 */
static uint64_t MASK_P;
static int BITS_PER;  /* bits per B value: ceil(log2(K)) */

/* Hash table */
#define HT_BITS 25   /* 32M slots = 512MB for p=43 */
#define HT_SIZE (1u << HT_BITS)
#define HT_MASK (HT_SIZE - 1)
#define EMPTY_HI 0xFFFFFFFFFFFFFFFFULL

typedef struct {
    uint64_t key_lo;
    uint64_t key_hi;
    uint32_t count;
    uint32_t _pad;
} ht_entry_t;

static ht_entry_t *ht;
static uint32_t ht_used = 0;

/* Marginal counts */
static uint64_t marginal[150][80];  /* marginal[j][b] for j < R, b < K */

static inline uint32_t ht_hash(uint64_t lo, uint64_t hi) {
    uint64_t h = lo;
    h ^= hi * 0x9e3779b97f4a7c15ULL;
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccdULL;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53ULL;
    h ^= h >> 33;
    return (uint32_t)(h & HT_MASK);
}

static inline void ht_insert(uint64_t lo, uint64_t hi) {
    uint32_t idx = ht_hash(lo, hi);
    while (1) {
        if (ht[idx].key_hi == EMPTY_HI && ht[idx].key_lo == EMPTY_HI) {
            ht[idx].key_lo = lo;
            ht[idx].key_hi = hi;
            ht[idx].count = 1;
            ht_used++;
            return;
        }
        if (ht[idx].key_lo == lo && ht[idx].key_hi == hi) {
            ht[idx].count++;
            return;
        }
        idx = (idx + 1) & HT_MASK;
    }
}

static inline uint64_t rot(uint64_t m, int d) {
    if (d == 0) return m;
    return ((m << d) | (m >> (P_VAL - d))) & MASK_P;
}

/* Gosper's hack for next combination */
static inline uint64_t next_comb(uint64_t v) {
    uint64_t t = v | (v - 1);
    return (t + 1) | (((~t & -(~t)) - 1) >> (__builtin_ctzll(v) + 1));
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <p>\n", argv[0]);
        return 1;
    }

    P_VAL = atoi(argv[1]);
    K_VAL = (P_VAL - 1) / 2;
    R_VAL = K_VAL;
    S_VAL = K_VAL * (K_VAL - 1) / 2;
    MASK_P = (1ULL << P_VAL) - 1;

    /* Determine bits per B value */
    BITS_PER = 0;
    int tmp = K_VAL;
    while (tmp > 0) { BITS_PER++; tmp >>= 1; }
    if (BITS_PER < 4) BITS_PER = 4;

    int total_bits = R_VAL * BITS_PER;

    printf("=== Min P/Q Verification at p=%d ===\n", P_VAL);
    printf("k=%d, R=%d, S=%d\n", K_VAL, R_VAL, S_VAL);
    printf("Bits per B value: %d, total bits: %d\n", BITS_PER, total_bits);

    if (P_VAL > 63) {
        fprintf(stderr, "ERROR: p=%d exceeds 63-bit mask limit\n", P_VAL);
        return 1;
    }
    if (total_bits > 128) {
        fprintf(stderr, "ERROR: B-profile needs %d bits > 128\n", total_bits);
        return 1;
    }

    /* Estimate total subsets */
    /* Use Stirling for C(p,k) */
    double log2_total = 0;
    for (int i = 0; i < K_VAL; i++)
        log2_total += log2((double)(P_VAL - i) / (i + 1));
    printf("log2(C(%d,%d)) = %.1f (approx %.3e subsets)\n",
           P_VAL, K_VAL, log2_total, pow(2.0, log2_total));

    printf("Hash table: %u slots (%lu MB)\n\n",
           HT_SIZE, (unsigned long)(HT_SIZE * sizeof(ht_entry_t) / 1000000));

    /* Allocate hash table */
    ht = (ht_entry_t *)malloc(HT_SIZE * sizeof(ht_entry_t));
    if (!ht) { fprintf(stderr, "OOM\n"); return 1; }
    for (uint32_t i = 0; i < HT_SIZE; i++) {
        ht[i].key_lo = EMPTY_HI;
        ht[i].key_hi = EMPTY_HI;
        ht[i].count = 0;
    }
    memset(marginal, 0, sizeof(marginal));

    time_t t0 = time(NULL);
    uint64_t total = 0;
    uint64_t d12 = (1ULL << K_VAL) - 1;

    while (1) {
        /* Compute B-profile and pack into 128-bit key */
        uint64_t key_lo = 0, key_hi = 0;
        int bit_pos = 0;
        for (int j = 0; j < R_VAL; j++) {
            int d = j + 1;
            uint64_t rd = rot(d12, d);
            int b = __builtin_popcountll(d12 & rd);

            /* Pack into key */
            if (bit_pos < 64) {
                key_lo |= ((uint64_t)b) << bit_pos;
                if (bit_pos + BITS_PER > 64) {
                    /* Straddles boundary */
                    key_hi |= ((uint64_t)b) >> (64 - bit_pos);
                }
            } else {
                key_hi |= ((uint64_t)b) << (bit_pos - 64);
            }
            bit_pos += BITS_PER;

            marginal[j][b]++;
        }

        ht_insert(key_lo, key_hi);
        total++;

        if ((total & 0xFFFFFF) == 0) {
            time_t now = time(NULL);
            long elapsed = now - t0;
            double total_est = pow(2.0, log2_total);
            double pct = total * 100.0 / total_est;
            double rate = elapsed > 0 ? total / (double)elapsed : 0;
            long eta = rate > 0 ? (long)((total_est - total) / rate) : 0;
            printf("  %12llu (%.2f%%) [%lds, ETA %lds] ht=%u (%.1f%%)\n",
                   (unsigned long long)total, pct, elapsed, eta,
                   ht_used, 100.0 * ht_used / HT_SIZE);
            fflush(stdout);

            /* Check hash table load */
            if (ht_used > HT_SIZE * 3 / 4) {
                fprintf(stderr, "ERROR: hash table >75%% full at %u entries. "
                        "Need larger table.\n", ht_used);
                break;
            }
        }

        d12 = next_comb(d12);
        if (d12 > MASK_P) break;
    }

    time_t t1 = time(NULL);
    printf("\nEnumeration complete: %llu subsets in %lds\n",
           (unsigned long long)total, t1 - t0);
    printf("Distinct B-profiles: %u\n", ht_used);
    printf("Hash table occupancy: %.1f%%\n\n", 100.0 * ht_used / HT_SIZE);

    /* Verify hyperplane */
    int hp_ok = 1;
    int checked = 0;
    for (uint32_t i = 0; i < HT_SIZE && checked < 100; i++) {
        if (ht[i].key_hi == EMPTY_HI && ht[i].key_lo == EMPTY_HI) continue;
        int sum = 0, bit_pos = 0;
        for (int j = 0; j < R_VAL; j++) {
            int b;
            if (bit_pos < 64) {
                b = (ht[i].key_lo >> bit_pos) & ((1 << BITS_PER) - 1);
            } else {
                b = (ht[i].key_hi >> (bit_pos - 64)) & ((1 << BITS_PER) - 1);
            }
            /* Handle straddling */
            if (bit_pos < 64 && bit_pos + BITS_PER > 64) {
                int lo_bits = 64 - bit_pos;
                b = (ht[i].key_lo >> bit_pos) & ((1 << lo_bits) - 1);
                b |= (ht[i].key_hi & ((1 << (BITS_PER - lo_bits)) - 1)) << lo_bits;
            }
            sum += b;
            bit_pos += BITS_PER;
        }
        if (sum != S_VAL) {
            printf("ERROR: profile sum=%d != %d\n", sum, S_VAL);
            hp_ok = 0;
            break;
        }
        checked++;
    }
    printf("Hyperplane check (%d profiles): %s\n", checked, hp_ok ? "PASS" : "FAIL");

    /* Check marginals identical */
    int marg_ok = 1;
    for (int j = 1; j < R_VAL; j++) {
        for (int b = 0; b < K_VAL; b++) {
            if (marginal[j][b] != marginal[0][b]) { marg_ok = 0; break; }
        }
        if (!marg_ok) break;
    }
    printf("Marginals identical: %s\n", marg_ok ? "YES" : "NO");

    /* Compute min P/Q */
    printf("\n=== P/Q RATIO ANALYSIS ===\n");
    double min_ratio = 1e30, max_ratio = 0;
    uint64_t below_1 = 0;
    double sum_log_ratio = 0;
    uint64_t n_profiles = 0;

    for (uint32_t i = 0; i < HT_SIZE; i++) {
        if (ht[i].key_hi == EMPTY_HI && ht[i].key_lo == EMPTY_HI) continue;

        double log_P = log((double)ht[i].count) - log((double)total);
        double log_Q = 0;
        int bit_pos = 0;
        for (int j = 0; j < R_VAL; j++) {
            int b;
            if (bit_pos < 64) {
                b = (ht[i].key_lo >> bit_pos) & ((1 << BITS_PER) - 1);
            } else {
                b = (ht[i].key_hi >> (bit_pos - 64)) & ((1 << BITS_PER) - 1);
            }
            if (bit_pos < 64 && bit_pos + BITS_PER > 64) {
                int lo_bits = 64 - bit_pos;
                b = (ht[i].key_lo >> bit_pos) & ((1 << lo_bits) - 1);
                b |= (ht[i].key_hi & ((1 << (BITS_PER - lo_bits)) - 1)) << lo_bits;
            }
            log_Q += log((double)marginal[j][b]) - log((double)total);
            bit_pos += BITS_PER;
        }

        double log_ratio = log_P - log_Q;
        double ratio = exp(log_ratio);

        if (ratio < min_ratio) min_ratio = ratio;
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
    }

    /* Write JSON results */
    char fname[256];
    snprintf(fname, sizeof(fname), "ramsey-book-graphs/min_pq_p%d_results.json", P_VAL);
    FILE *fout = fopen(fname, "w");
    if (fout) {
        fprintf(fout, "{\n");
        fprintf(fout, "  \"p\": %d,\n", P_VAL);
        fprintf(fout, "  \"k\": %d,\n", K_VAL);
        fprintf(fout, "  \"total_subsets\": %llu,\n", (unsigned long long)total);
        fprintf(fout, "  \"distinct_profiles\": %u,\n", ht_used);
        fprintf(fout, "  \"min_PQ\": %.10f,\n", min_ratio);
        fprintf(fout, "  \"max_PQ\": %.6e,\n", max_ratio);
        fprintf(fout, "  \"below_1\": %llu,\n", (unsigned long long)below_1);
        fprintf(fout, "  \"min_PQ_gt_1\": %s,\n", min_ratio > 1.0 ? "true" : "false");
        fprintf(fout, "  \"log2_min_PQ\": %.6f\n", log(min_ratio) / log(2.0));
        fprintf(fout, "}\n");
        fclose(fout);
        printf("\nResults written to %s\n", fname);
    }

    free(ht);
    return 0;
}
