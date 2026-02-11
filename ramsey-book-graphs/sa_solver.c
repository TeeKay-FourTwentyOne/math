/*
 * High-performance SA solver for R(B_{n-1}, B_n) = 4n-1.
 *
 * Uses O(m) incremental delta updates per move.
 * Target throughput: >1M iterations/sec for m=63.
 *
 * Usage: ./sa_solver <n> [max_seeds] [max_iter]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_M 200

/* Global parameters */
static int N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D12_SIZE;

/* State arrays (indicators, 0 or 1) */
static int D11[MAX_M];
static int D12[MAX_M];
static int D12T[MAX_M];

/* Autocorrelation: delta_XX[d] = |{a in X : (a-d) mod m in X}| */
static int delta_11[MAX_M];
static int delta_12[MAX_M];
static int delta_12T[MAX_M];

/* V1V2 verification (only checked when V1V1/V2V2 cost = 0) */

/* Symmetric pairs for D11 */
static int pairs[MAX_M][2];
static int num_pairs;
static int pair_in_d11[MAX_M]; /* 1 if pair i is currently in D11 */

/* Lists for O(1) random selection in D12 */
static int d12_mem[MAX_M], d12_nmem[MAX_M];
static int d12_cnt, d12_ncnt;
/* Index arrays for O(1) list updates */
static int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];

/* Pair lists for O(1) random selection in D11 */
static int d11_in[MAX_M], d11_out[MAX_M];
static int d11_in_cnt, d11_out_cnt;
static int d11_in_idx[MAX_M], d11_out_idx[MAX_M];

/* Best state */
static int best_D11[MAX_M], best_D12[MAX_M];
static int best_cost;

/* Precomputed modular arithmetic tables */
static int mod_add[MAX_M][MAX_M]; /* (a+b) % M */
static int mod_sub[MAX_M][MAX_M]; /* (a-b+M) % M */
static int mod_neg[MAX_M];        /* (M-a) % M */

/* Random number generator (xoshiro256**) */
static unsigned long long rng_s[4];

static inline unsigned long long rotl(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(void) {
    const unsigned long long result = rotl(rng_s[1] * 5, 7) * 9;
    const unsigned long long t = rng_s[1] << 17;
    rng_s[2] ^= rng_s[0];
    rng_s[3] ^= rng_s[1];
    rng_s[1] ^= rng_s[2];
    rng_s[0] ^= rng_s[3];
    rng_s[2] ^= t;
    rng_s[3] = rotl(rng_s[3], 45);
    return result;
}

static void rng_seed(unsigned long long seed) {
    rng_s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    rng_s[1] = seed ^ 0x6C62272E07BB0142ULL;
    rng_s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    rng_s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next();
}

static inline int rng_int(int n) {
    return (int)(rng_next() % (unsigned long long)n);
}

static inline double rng_double(void) {
    return (rng_next() >> 11) * 0x1.0p-53;
}

/* Initialize modular arithmetic tables */
static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add[a][b] = (a + b) % M;
            mod_sub[a][b] = (a - b + M) % M;
        }
    }
}

/* Compute full autocorrelation from scratch */
static void compute_delta_full(const int *ind, int *delta) {
    for (int d = 0; d < M; d++) {
        int count = 0;
        for (int a = 0; a < M; a++) {
            if (ind[a] && ind[mod_sub[a][d]])
                count++;
        }
        delta[d] = count;
    }
}

/*
 * Incremental delta update when swapping rem -> add_el in set S.
 * S indicator must already reflect the new state.
 * Updates delta[] in-place. O(m) per call.
 */
static inline void update_delta_inc(int *delta, const int *S, int rem, int add_el) {
    /* delta[0] = |S| doesn't change (swap preserves size).
     *
     * For d != 0, the change is:
     *   new_delta[d] - old_delta[d] =
     *     S_new[add-d] + S_new[add+d] - S_old[rem-d] - S_old[rem+d]
     *
     * Since S_old(x) = S_new(x) + [x==rem] - [x==add_el], we get:
     *   change = S_new[add-d] + S_new[add+d] - S_new[rem-d] - S_new[rem+d]
     *          + [(rem-d)==add_el] + [(rem+d)==add_el]
     */
    for (int d = 1; d < M; d++) {
        int change = S[mod_sub[add_el][d]] + S[mod_add[add_el][d]]
                   - S[mod_sub[rem][d]] - S[mod_add[rem][d]];
        if (mod_sub[rem][d] == add_el) change++;
        if (mod_add[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

/* Revert an incremental delta update (undo swap of rem -> add_el).
 * S indicator must reflect the REVERTED state (rem back in, add_el out). */
static inline void revert_delta_inc(int *delta, const int *S, int rem, int add_el) {
    /* This is equivalent to update_delta_inc with rem and add_el swapped */
    update_delta_inc(delta, S, add_el, rem);
}

/* Compute cost from current delta arrays and indicators (V1V1 + V2V2 only) */
static int compute_cost(void) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += D11[i];
        d12_size += D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int cost = 0;

    for (int d = 1; d < M; d++) {
        int cv11 = delta_11[d] + delta_12[d];
        if (D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        /* V2V2 via complement formula */
        int d22d = delta_11[d] + (M - 2 - 2 * d11_size) + 2 * D11[d];
        int cv22 = d22d + delta_12T[d];
        if (!D11[d]) {
            int excess = cv22 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }

    return cost;
}

/* Full V1V2 verification (O(m^2), only called when V1V1/V2V2 cost=0) */
static int verify_v1v2(void) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += D11[i];
        d12_size += D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[MAX_M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        /* Sigma(D11, D12, d) = |{a in D11 : (d-a) mod m in D12}| */
        int sigma = 0;
        for (int a = 0; a < M; a++) {
            if (D11[a] && D12[mod_sub[d][a]])
                sigma++;
        }
        /* Delta(D12, D22, d) = |{a in D12 : (a-d) mod m in D22}| */
        int delta = 0;
        for (int a = 0; a < M; a++) {
            if (D12[a] && D22[mod_sub[a][d]])
                delta++;
        }
        int common = sigma + delta;

        if (D12[d]) {
            /* Red cross-edge */
            int excess = common - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            /* Blue cross-edge */
            int bc = (N_TOTAL - 2) - d1 - d2 + common;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* Build D12 member/nonmember lists with index tracking */
static void rebuild_d12_lists(void) {
    d12_cnt = 0;
    d12_ncnt = 0;
    memset(d12_mem_idx, -1, sizeof(d12_mem_idx));
    memset(d12_nmem_idx, -1, sizeof(d12_nmem_idx));
    for (int i = 0; i < M; i++) {
        if (D12[i]) {
            d12_mem_idx[i] = d12_cnt;
            d12_mem[d12_cnt++] = i;
        } else {
            d12_nmem_idx[i] = d12_ncnt;
            d12_nmem[d12_ncnt++] = i;
        }
    }
}

/* Build D11 pair in/out lists with index tracking */
static void rebuild_d11_lists(void) {
    d11_in_cnt = 0;
    d11_out_cnt = 0;
    for (int i = 0; i < num_pairs; i++) {
        if (pair_in_d11[i]) {
            d11_in_idx[i] = d11_in_cnt;
            d11_in[d11_in_cnt++] = i;
        } else {
            d11_out_idx[i] = d11_out_cnt;
            d11_out[d11_out_cnt++] = i;
        }
    }
}

/* Swap elements in D12 list tracking */
static void d12_swap_to_mem(int elem) {
    /* Move elem from nonmember to member list */
    int ni = d12_nmem_idx[elem];
    /* Replace with last element in nonmember list */
    int last = d12_nmem[d12_ncnt - 1];
    d12_nmem[ni] = last;
    d12_nmem_idx[last] = ni;
    d12_ncnt--;
    /* Add to member list */
    d12_mem_idx[elem] = d12_cnt;
    d12_mem[d12_cnt++] = elem;
    d12_nmem_idx[elem] = -1;
}

static void d12_swap_to_nmem(int elem) {
    int mi = d12_mem_idx[elem];
    int last = d12_mem[d12_cnt - 1];
    d12_mem[mi] = last;
    d12_mem_idx[last] = mi;
    d12_cnt--;
    d12_nmem_idx[elem] = d12_ncnt;
    d12_nmem[d12_ncnt++] = elem;
    d12_mem_idx[elem] = -1;
}

/* Build symmetric pairs */
static void build_pairs(void) {
    num_pairs = 0;
    for (int x = 1; x < M; x++) {
        int neg_x = mod_neg[x];
        if (x <= neg_x) {
            pairs[num_pairs][0] = x;
            pairs[num_pairs][1] = neg_x;
            num_pairs++;
        }
    }
}

static int solve(int d11_size, int seed, int max_iter) {
    rng_seed((unsigned long long)seed * 1000003ULL + 42);

    /* Initialize D11 randomly */
    memset(D11, 0, sizeof(int) * M);
    int num_selected = d11_size / 2;

    int perm[MAX_M];
    for (int i = 0; i < num_pairs; i++) perm[i] = i;
    for (int i = num_pairs - 1; i > 0; i--) {
        int j = rng_int(i + 1);
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    for (int i = 0; i < num_selected; i++) {
        D11[pairs[perm[i]][0]] = 1;
        D11[pairs[perm[i]][1]] = 1;
    }
    memset(pair_in_d11, 0, sizeof(int) * num_pairs);
    for (int i = 0; i < num_selected; i++)
        pair_in_d11[perm[i]] = 1;
    rebuild_d11_lists();

    /* Initialize D12 randomly - 0 is NOT fixed in D12.
     * Select D12_SIZE elements from {0, 1, ..., M-1} randomly. */
    memset(D12, 0, sizeof(int) * M);
    memset(D12T, 0, sizeof(int) * M);

    {
        int elems[MAX_M];
        for (int i = 0; i < M; i++) elems[i] = i;
        for (int i = M - 1; i > 0; i--) {
            int j = rng_int(i + 1);
            int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
        }
        for (int i = 0; i < D12_SIZE; i++) {
            D12[elems[i]] = 1;
            D12T[mod_neg[elems[i]]] = 1;
        }
    }
    rebuild_d12_lists();

    /* Compute initial deltas from scratch */
    compute_delta_full(D11, delta_11);
    compute_delta_full(D12, delta_12);
    compute_delta_full(D12T, delta_12T);

    int current_cost = compute_cost();
    best_cost = current_cost;
    memcpy(best_D11, D11, sizeof(int) * M);
    memcpy(best_D12, D12, sizeof(int) * M);

    /* Very slow cooling with no reheat. Cool from T_start to T_min over max_iter.
     * alpha = exp(ln(T_min/T_start) / max_iter) */
    double T_start = 15.0;
    double T_min = 0.0001;
    double T = T_start;
    double alpha = exp(log(T_min / T_start) / (double)max_iter);
    int last_imp = 0;
    int reheat_count = 0;

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    /* Precompute log table for acceptance (optional speedup) */

    for (int it = 0; it < max_iter; it++) {
        int r = rng_int(100);

        if (r < 22) {
            /* D11 pair swap */
            if (d11_in_cnt == 0 || d11_out_cnt == 0) continue;

            int ri = rng_int(d11_in_cnt);
            int oi = rng_int(d11_out_cnt);
            int rp = d11_in[ri]; /* pair index to remove */
            int ap = d11_out[oi]; /* pair index to add */

            int r0 = pairs[rp][0], r1 = pairs[rp][1];
            int a0 = pairs[ap][0], a1 = pairs[ap][1];

            /* Apply move to D11 */
            D11[r0] = 0; D11[r1] = 0;
            D11[a0] = 1; D11[a1] = 1;

            /* For D11, we swap TWO elements at once (a symmetric pair).
             * We need to update delta_11 incrementally.
             * First remove r0, then remove r1, then add a0, then add a1.
             * Actually, we can decompose into two single-element swaps... but that's
             * not quite right since it's removing 2 and adding 2.
             *
             * For a 2-element removal + 2-element addition, it's easier to just
             * recompute delta_11 from scratch. For m=63, this is 63*63 = 3969 ops.
             */
            compute_delta_full(D11, delta_11);
            int new_cost = compute_cost();
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                pair_in_d11[rp] = 0;
                pair_in_d11[ap] = 1;
                rebuild_d11_lists();
            } else {
                /* Revert */
                D11[r0] = 1; D11[r1] = 1;
                D11[a0] = 0; D11[a1] = 0;
                compute_delta_full(D11, delta_11);
            }
        } else if (r < 82) {
            /* D12 single swap */
            if (d12_cnt == 0 || d12_ncnt == 0) continue;

            int ri = rng_int(d12_cnt);
            int ai = rng_int(d12_ncnt);
            int rem = d12_mem[ri];
            int add_el = d12_nmem[ai];

            /* Apply swap */
            D12[rem] = 0;
            D12[add_el] = 1;
            int rem_t = mod_neg[rem];
            int add_t = mod_neg[add_el];
            D12T[rem_t] = 0;
            D12T[add_t] = 1;

            /* Incremental delta update */
            update_delta_inc(delta_12, D12, rem, add_el);
            update_delta_inc(delta_12T, D12T, rem_t, add_t);

            int new_cost = compute_cost();
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(rem);
                d12_swap_to_mem(add_el);
            } else {
                /* Revert */
                D12[rem] = 1;
                D12[add_el] = 0;
                D12T[rem_t] = 1;
                D12T[add_t] = 0;
                revert_delta_inc(delta_12, D12, rem, add_el);
                revert_delta_inc(delta_12T, D12T, rem_t, add_t);
            }
        } else {
            /* D12 double swap */
            if (d12_cnt < 2 || d12_ncnt < 2) continue;

            int ri1 = rng_int(d12_cnt);
            int ri2 = rng_int(d12_cnt - 1);
            if (ri2 >= ri1) ri2++;
            int ai1 = rng_int(d12_ncnt);
            int ai2 = rng_int(d12_ncnt - 1);
            if (ai2 >= ai1) ai2++;

            int rem1 = d12_mem[ri1], rem2 = d12_mem[ri2];
            int add1 = d12_nmem[ai1], add2 = d12_nmem[ai2];

            /* For double swap, recompute from scratch (easier, still fast) */
            int old_d12[MAX_M], old_d12T[MAX_M];
            memcpy(old_d12, delta_12, sizeof(int) * M);
            memcpy(old_d12T, delta_12T, sizeof(int) * M);

            D12[rem1] = 0; D12[rem2] = 0;
            D12[add1] = 1; D12[add2] = 1;
            D12T[mod_neg[rem1]] = 0; D12T[mod_neg[rem2]] = 0;
            D12T[mod_neg[add1]] = 1; D12T[mod_neg[add2]] = 1;

            compute_delta_full(D12, delta_12);
            compute_delta_full(D12T, delta_12T);

            int new_cost = compute_cost();
            int dc = new_cost - current_cost;

            if (dc <= 0 || rng_double() < exp(-(double)dc / fmax(T, T_min))) {
                current_cost = new_cost;
                d12_swap_to_nmem(rem1);
                d12_swap_to_nmem(rem2);
                d12_swap_to_mem(add1);
                d12_swap_to_mem(add2);
            } else {
                D12[rem1] = 1; D12[rem2] = 1;
                D12[add1] = 0; D12[add2] = 0;
                D12T[mod_neg[rem1]] = 1; D12T[mod_neg[rem2]] = 1;
                D12T[mod_neg[add1]] = 0; D12T[mod_neg[add2]] = 0;
                memcpy(delta_12, old_d12, sizeof(int) * M);
                memcpy(delta_12T, old_d12T, sizeof(int) * M);
            }
        }

        if (current_cost < best_cost) {
            best_cost = current_cost;
            memcpy(best_D11, D11, sizeof(int) * M);
            memcpy(best_D12, D12, sizeof(int) * M);
            last_imp = it;
            if (best_cost == 0) {
                /* Check V1V2 constraints too */
                int v12_cost = verify_v1v2();
                if (v12_cost == 0) {
                    break; /* Fully valid! */
                } else {
                    /* V1V1/V2V2 OK but V1V2 fails - include V1V2 in cost */
                    best_cost = v12_cost;
                }
            }
        }

        T *= alpha;

        if (it % 2000000 == 0 && it > 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                           (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
            printf("    %dM: cost=%d best=%d T=%.5f rh=%d (%.0f it/s)\n",
                   it / 1000000, current_cost, best_cost, T, reheat_count,
                   (double)it / elapsed);
            fflush(stdout);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                   (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("  seed=%d: best=%d (%.1fs)\n", seed, best_cost, elapsed);

    /* Dump diagnostic info for near-solutions */
    if (best_cost > 0 && best_cost <= 8) {
        /* Restore best state */
        memcpy(D11, best_D11, sizeof(int) * M);
        memcpy(D12, best_D12, sizeof(int) * M);
        memset(D12T, 0, sizeof(int) * M);
        for (int i = 0; i < M; i++)
            if (D12[i]) D12T[mod_neg[i]] = 1;
        compute_delta_full(D11, delta_11);
        compute_delta_full(D12, delta_12);
        compute_delta_full(D12T, delta_12T);

        int d11_sz = 0, d12_sz = 0;
        for (int i = 0; i < M; i++) { d11_sz += D11[i]; d12_sz += D12[i]; }
        int d22_sz = M - 1 - d11_sz;
        int d1 = d11_sz + d12_sz, d2 = d22_sz + d12_sz;

        printf("  VIOLATIONS (V1V1+V2V2):\n");
        for (int d = 1; d < M; d++) {
            int cv11 = delta_11[d] + delta_12[d];
            if (D11[d]) {
                int excess = cv11 - RED_THRESH;
                if (excess > 0) printf("    V1V1 red d=%d: cv=%d thresh=%d excess=%d\n", d, cv11, RED_THRESH, excess);
            } else {
                int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
                int excess = bc - BLUE_THRESH;
                if (excess > 0) printf("    V1V1 blue d=%d: bc=%d thresh=%d excess=%d\n", d, bc, BLUE_THRESH, excess);
            }
            int d22d = delta_11[d] + (M - 2 - 2 * d11_sz) + 2 * D11[d];
            int cv22 = d22d + delta_12T[d];
            if (!D11[d]) {
                int excess = cv22 - RED_THRESH;
                if (excess > 0) printf("    V2V2 red d=%d: cv=%d thresh=%d excess=%d\n", d, cv22, RED_THRESH, excess);
            } else {
                int bc = (N_TOTAL - 2) - 2 * d2 + cv22;
                int excess = bc - BLUE_THRESH;
                if (excess > 0) printf("    V2V2 blue d=%d: bc=%d thresh=%d excess=%d\n", d, bc, BLUE_THRESH, excess);
            }
        }

        /* Also check V1V2 on this state */
        int v12 = verify_v1v2();
        printf("  V1V2 cost on best state: %d\n", v12);

        /* Dump D11 and D12 as JSON for analysis */
        if (best_cost == 4) {
            char fname[256];
            snprintf(fname, sizeof(fname), "/tmp/near_n%d_seed%d.json", N_PARAM, seed);
            FILE *fp = fopen(fname, "w");
            if (fp) {
                fprintf(fp, "{\"n\": %d, \"m\": %d, \"cost\": %d,\n", N_PARAM, M, best_cost);
                fprintf(fp, " \"D11\": [");
                for (int i = 0, first = 1; i < M; i++)
                    if (best_D11[i]) { if (!first) fprintf(fp, ", "); fprintf(fp, "%d", i); first = 0; }
                fprintf(fp, "],\n \"D12\": [");
                for (int i = 0, first = 1; i < M; i++)
                    if (best_D12[i]) { if (!first) fprintf(fp, ", "); fprintf(fp, "%d", i); first = 0; }
                fprintf(fp, "]}\n");
                fclose(fp);
                printf("  Dumped near-solution to %s\n", fname);
            }
        }
    }

    fflush(stdout);

    return best_cost;
}


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [max_seeds] [max_iter] [fixed_d11_size]\n", argv[0]);
        return 1;
    }

    N_PARAM = atoi(argv[1]);
    int max_seeds = argc > 2 ? atoi(argv[2]) : 50;
    int max_iter = argc > 3 ? atoi(argv[3]) : 50000000;
    int fixed_d11 = argc > 4 ? atoi(argv[4]) : -1;

    M = 2 * N_PARAM - 1;
    N_TOTAL = 2 * M;
    RED_THRESH = N_PARAM - 2;
    BLUE_THRESH = N_PARAM - 1;
    D12_SIZE = N_PARAM - 1;

    printf("n=%d, m=%d, N=%d\n", N_PARAM, M, N_TOTAL);
    printf("red_thresh=%d, blue_thresh=%d, |D12|=%d\n", RED_THRESH, BLUE_THRESH, D12_SIZE);
    printf("Seeds: %d, Iter: %d\n", max_seeds, max_iter);
    fflush(stdout);

    init_mod_tables();
    build_pairs();
    printf("Symmetric pairs: %d\n", num_pairs);
    fflush(stdout);

    /* Determine D11 sizes to try */
    int half = (M - 1) / 2;
    int d11_sizes[10];
    int num_sizes = 0;

    if (fixed_d11 > 0) {
        d11_sizes[0] = fixed_d11;
        num_sizes = 1;
    } else {
        for (int offset = 0; offset <= 3; offset++) {
            int candidates[2] = {half - offset, half + offset};
            for (int c = 0; c < 2; c++) {
                int s = candidates[c];
                if (s > 0 && s < M - 1 && s % 2 == 0 && (s / 2) <= num_pairs) {
                    int dup = 0;
                    for (int j = 0; j < num_sizes; j++)
                        if (d11_sizes[j] == s) { dup = 1; break; }
                    if (!dup && num_sizes < 10) d11_sizes[num_sizes++] = s;
                }
            }
        }

        /* Sort by distance from half */
        for (int i = 0; i < num_sizes; i++) {
            for (int j = i + 1; j < num_sizes; j++) {
                if (abs(d11_sizes[j] - half) < abs(d11_sizes[i] - half)) {
                    int tmp = d11_sizes[i];
                    d11_sizes[i] = d11_sizes[j];
                    d11_sizes[j] = tmp;
                }
            }
        }
    }

    printf("D11 sizes to try:");
    for (int i = 0; i < num_sizes; i++) printf(" %d", d11_sizes[i]);
    printf("\n");
    fflush(stdout);

    int global_best = 999;

    for (int si = 0; si < num_sizes; si++) {
        int d11_size = d11_sizes[si];
        printf("\n--- |D11| = %d ---\n", d11_size);
        fflush(stdout);

        for (int seed = 0; seed < max_seeds; seed++) {
            int result = solve(d11_size, seed, max_iter);
            if (result < global_best) global_best = result;

            if (result == 0) {
                printf("\n*** SOLUTION FOUND! ***\n");
                printf("n=%d, m=%d\n", N_PARAM, M);
                printf("D11 = [");
                for (int i = 0, first = 1; i < M; i++)
                    if (best_D11[i]) { if (!first) printf(", "); printf("%d", i); first = 0; }
                printf("]\n");
                printf("D12 = [");
                for (int i = 0, first = 1; i < M; i++)
                    if (best_D12[i]) { if (!first) printf(", "); printf("%d", i); first = 0; }
                printf("]\n");
                printf("|D11|=%d, |D12|=%d\n", d11_size, D12_SIZE);
                fflush(stdout);
                return 0;
            }
        }
    }

    printf("\nGlobal best = %d\n", global_best);
    return 1;
}
