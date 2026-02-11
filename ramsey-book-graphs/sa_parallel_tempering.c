/*
 * Parallel Tempering (Replica Exchange) SA solver for R(B_{n-1}, B_n) = 4n-1.
 *
 * Runs R replicas at different temperatures simultaneously.
 * Periodically attempts to swap adjacent replicas using Metropolis criterion.
 * This helps escape local minima that trap standard SA at cost=4.
 *
 * Usage: ./sa_pt <n> [max_seeds] [max_iter] [fixed_d11_size]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_M 200
#define NUM_REPLICAS 16
#define SWAP_INTERVAL 500  /* Attempt replica swaps every N iterations */

/* Global parameters */
static int N_PARAM, M, N_TOTAL;
static int RED_THRESH, BLUE_THRESH;
static int D12_SIZE;

/* Precomputed modular arithmetic tables */
static int mod_add[MAX_M][MAX_M];
static int mod_sub[MAX_M][MAX_M];
static int mod_neg[MAX_M];

/* Symmetric pairs for D11 */
static int pairs[MAX_M][2];
static int num_pairs;

/* Per-replica state */
typedef struct {
    int D11[MAX_M];
    int D12[MAX_M];
    int D12T[MAX_M];
    int delta_11[MAX_M];
    int delta_12[MAX_M];
    int delta_12T[MAX_M];
    int pair_in_d11[MAX_M];

    /* D12 member/nonmember lists */
    int d12_mem[MAX_M], d12_nmem[MAX_M];
    int d12_cnt, d12_ncnt;
    int d12_mem_idx[MAX_M], d12_nmem_idx[MAX_M];

    /* D11 pair in/out lists */
    int d11_in[MAX_M], d11_out[MAX_M];
    int d11_in_cnt, d11_out_cnt;
    int d11_in_idx[MAX_M], d11_out_idx[MAX_M];

    double temperature;
    int current_cost;

    /* RNG state (per-replica) */
    unsigned long long rng_s[4];
} Replica;

static Replica replicas[NUM_REPLICAS];

/* Best solution found globally */
static int best_D11[MAX_M], best_D12[MAX_M];
static int global_best_cost;

/* --- RNG (per-replica) --- */

static inline unsigned long long rotl(const unsigned long long x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline unsigned long long rng_next(Replica *r) {
    const unsigned long long result = rotl(r->rng_s[1] * 5, 7) * 9;
    const unsigned long long t = r->rng_s[1] << 17;
    r->rng_s[2] ^= r->rng_s[0];
    r->rng_s[3] ^= r->rng_s[1];
    r->rng_s[1] ^= r->rng_s[2];
    r->rng_s[0] ^= r->rng_s[3];
    r->rng_s[2] ^= t;
    r->rng_s[3] = rotl(r->rng_s[3], 45);
    return result;
}

static void rng_seed(Replica *r, unsigned long long seed) {
    r->rng_s[0] = seed ^ 0x9E3779B97F4A7C15ULL;
    r->rng_s[1] = seed ^ 0x6C62272E07BB0142ULL;
    r->rng_s[2] = seed ^ 0xBF58476D1CE4E5B9ULL;
    r->rng_s[3] = seed ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) rng_next(r);
}

static inline int rng_int(Replica *r, int n) {
    return (int)(rng_next(r) % (unsigned long long)n);
}

static inline double rng_double(Replica *r) {
    return (rng_next(r) >> 11) * 0x1.0p-53;
}

/* --- Modular arithmetic --- */

static void init_mod_tables(void) {
    for (int a = 0; a < M; a++) {
        mod_neg[a] = (M - a) % M;
        for (int b = 0; b < M; b++) {
            mod_add[a][b] = (a + b) % M;
            mod_sub[a][b] = (a - b + M) % M;
        }
    }
}

/* --- Delta computation --- */

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

static inline void update_delta_inc(int *delta, const int *S, int rem, int add_el) {
    for (int d = 1; d < M; d++) {
        int change = S[mod_sub[add_el][d]] + S[mod_add[add_el][d]]
                   - S[mod_sub[rem][d]] - S[mod_add[rem][d]];
        if (mod_sub[rem][d] == add_el) change++;
        if (mod_add[rem][d] == add_el) change++;
        delta[d] += change;
    }
}

static inline void revert_delta_inc(int *delta, const int *S, int rem, int add_el) {
    update_delta_inc(delta, S, add_el, rem);
}

/* --- Cost function (V1V1 + V2V2) --- */

static int compute_cost_r(Replica *r) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += r->D11[i];
        d12_size += r->D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int cost = 0;
    for (int d = 1; d < M; d++) {
        int cv11 = r->delta_11[d] + r->delta_12[d];
        if (r->D11[d]) {
            int excess = cv11 - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - 2 * d1 + cv11;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }

        int d22d = r->delta_11[d] + (M - 2 - 2 * d11_size) + 2 * r->D11[d];
        int cv22 = d22d + r->delta_12T[d];
        if (!r->D11[d]) {
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

/* --- V1V2 verification --- */

static int verify_v1v2(Replica *r) {
    int d11_size = 0, d12_size = 0;
    for (int i = 0; i < M; i++) {
        d11_size += r->D11[i];
        d12_size += r->D12[i];
    }
    int d22_size = M - 1 - d11_size;
    int d1 = d11_size + d12_size;
    int d2 = d22_size + d12_size;

    int D22[MAX_M];
    for (int i = 0; i < M; i++)
        D22[i] = (i > 0 && !r->D11[i]) ? 1 : 0;

    int cost = 0;
    for (int d = 0; d < M; d++) {
        int sigma = 0;
        for (int a = 0; a < M; a++) {
            if (r->D11[a] && r->D12[mod_sub[d][a]])
                sigma++;
        }
        int delta = 0;
        for (int a = 0; a < M; a++) {
            if (r->D12[a] && D22[mod_sub[a][d]])
                delta++;
        }
        int common = sigma + delta;

        if (r->D12[d]) {
            int excess = common - RED_THRESH;
            if (excess > 0) cost += excess;
        } else {
            int bc = (N_TOTAL - 2) - d1 - d2 + common;
            int excess = bc - BLUE_THRESH;
            if (excess > 0) cost += excess;
        }
    }
    return cost;
}

/* --- List management --- */

static void rebuild_d12_lists(Replica *r) {
    r->d12_cnt = 0;
    r->d12_ncnt = 0;
    memset(r->d12_mem_idx, -1, sizeof(r->d12_mem_idx));
    memset(r->d12_nmem_idx, -1, sizeof(r->d12_nmem_idx));
    for (int i = 0; i < M; i++) {
        if (r->D12[i]) {
            r->d12_mem_idx[i] = r->d12_cnt;
            r->d12_mem[r->d12_cnt++] = i;
        } else {
            r->d12_nmem_idx[i] = r->d12_ncnt;
            r->d12_nmem[r->d12_ncnt++] = i;
        }
    }
}

static void rebuild_d11_lists(Replica *r) {
    r->d11_in_cnt = 0;
    r->d11_out_cnt = 0;
    for (int i = 0; i < num_pairs; i++) {
        if (r->pair_in_d11[i]) {
            r->d11_in_idx[i] = r->d11_in_cnt;
            r->d11_in[r->d11_in_cnt++] = i;
        } else {
            r->d11_out_idx[i] = r->d11_out_cnt;
            r->d11_out[r->d11_out_cnt++] = i;
        }
    }
}

static void d12_swap_to_mem(Replica *r, int elem) {
    int ni = r->d12_nmem_idx[elem];
    int last = r->d12_nmem[r->d12_ncnt - 1];
    r->d12_nmem[ni] = last;
    r->d12_nmem_idx[last] = ni;
    r->d12_ncnt--;
    r->d12_mem_idx[elem] = r->d12_cnt;
    r->d12_mem[r->d12_cnt++] = elem;
    r->d12_nmem_idx[elem] = -1;
}

static void d12_swap_to_nmem(Replica *r, int elem) {
    int mi = r->d12_mem_idx[elem];
    int last = r->d12_mem[r->d12_cnt - 1];
    r->d12_mem[mi] = last;
    r->d12_mem_idx[last] = mi;
    r->d12_cnt--;
    r->d12_nmem_idx[elem] = r->d12_ncnt;
    r->d12_nmem[r->d12_ncnt++] = elem;
    r->d12_mem_idx[elem] = -1;
}

/* --- Pair management --- */

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

/* --- Initialize a replica --- */

static void init_replica(Replica *r, int d11_size, int seed, int replica_idx) {
    rng_seed(r, (unsigned long long)seed * 1000003ULL + 42 + replica_idx * 999983ULL);

    /* Initialize D11 randomly */
    memset(r->D11, 0, sizeof(int) * M);
    int num_selected = d11_size / 2;

    int perm[MAX_M];
    for (int i = 0; i < num_pairs; i++) perm[i] = i;
    for (int i = num_pairs - 1; i > 0; i--) {
        int j = rng_int(r, i + 1);
        int tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp;
    }
    for (int i = 0; i < num_selected; i++) {
        r->D11[pairs[perm[i]][0]] = 1;
        r->D11[pairs[perm[i]][1]] = 1;
    }
    memset(r->pair_in_d11, 0, sizeof(int) * num_pairs);
    for (int i = 0; i < num_selected; i++)
        r->pair_in_d11[perm[i]] = 1;
    rebuild_d11_lists(r);

    /* Initialize D12 randomly */
    memset(r->D12, 0, sizeof(int) * M);
    memset(r->D12T, 0, sizeof(int) * M);

    {
        int elems[MAX_M];
        for (int i = 0; i < M; i++) elems[i] = i;
        for (int i = M - 1; i > 0; i--) {
            int j = rng_int(r, i + 1);
            int tmp = elems[i]; elems[i] = elems[j]; elems[j] = tmp;
        }
        for (int i = 0; i < D12_SIZE; i++) {
            r->D12[elems[i]] = 1;
            r->D12T[mod_neg[elems[i]]] = 1;
        }
    }
    rebuild_d12_lists(r);

    compute_delta_full(r->D11, r->delta_11);
    compute_delta_full(r->D12, r->delta_12);
    compute_delta_full(r->D12T, r->delta_12T);

    r->current_cost = compute_cost_r(r);
}

/* --- One SA step on a replica --- */

static void sa_step(Replica *r) {
    int rr = rng_int(r, 100);

    if (rr < 22) {
        /* D11 pair swap */
        if (r->d11_in_cnt == 0 || r->d11_out_cnt == 0) return;

        int ri = rng_int(r, r->d11_in_cnt);
        int oi = rng_int(r, r->d11_out_cnt);
        int rp = r->d11_in[ri];
        int ap = r->d11_out[oi];

        int r0 = pairs[rp][0], r1 = pairs[rp][1];
        int a0 = pairs[ap][0], a1 = pairs[ap][1];

        r->D11[r0] = 0; r->D11[r1] = 0;
        r->D11[a0] = 1; r->D11[a1] = 1;

        compute_delta_full(r->D11, r->delta_11);
        int new_cost = compute_cost_r(r);
        int dc = new_cost - r->current_cost;

        if (dc <= 0 || rng_double(r) < exp(-(double)dc / fmax(r->temperature, 1e-10))) {
            r->current_cost = new_cost;
            r->pair_in_d11[rp] = 0;
            r->pair_in_d11[ap] = 1;
            rebuild_d11_lists(r);
        } else {
            r->D11[r0] = 1; r->D11[r1] = 1;
            r->D11[a0] = 0; r->D11[a1] = 0;
            compute_delta_full(r->D11, r->delta_11);
        }
    } else if (rr < 82) {
        /* D12 single swap */
        if (r->d12_cnt == 0 || r->d12_ncnt == 0) return;

        int ri = rng_int(r, r->d12_cnt);
        int ai = rng_int(r, r->d12_ncnt);
        int rem = r->d12_mem[ri];
        int add_el = r->d12_nmem[ai];

        r->D12[rem] = 0;
        r->D12[add_el] = 1;
        int rem_t = mod_neg[rem];
        int add_t = mod_neg[add_el];
        r->D12T[rem_t] = 0;
        r->D12T[add_t] = 1;

        update_delta_inc(r->delta_12, r->D12, rem, add_el);
        update_delta_inc(r->delta_12T, r->D12T, rem_t, add_t);

        int new_cost = compute_cost_r(r);
        int dc = new_cost - r->current_cost;

        if (dc <= 0 || rng_double(r) < exp(-(double)dc / fmax(r->temperature, 1e-10))) {
            r->current_cost = new_cost;
            d12_swap_to_nmem(r, rem);
            d12_swap_to_mem(r, add_el);
        } else {
            r->D12[rem] = 1;
            r->D12[add_el] = 0;
            r->D12T[rem_t] = 1;
            r->D12T[add_t] = 0;
            revert_delta_inc(r->delta_12, r->D12, rem, add_el);
            revert_delta_inc(r->delta_12T, r->D12T, rem_t, add_t);
        }
    } else {
        /* D12 double swap */
        if (r->d12_cnt < 2 || r->d12_ncnt < 2) return;

        int ri1 = rng_int(r, r->d12_cnt);
        int ri2 = rng_int(r, r->d12_cnt - 1);
        if (ri2 >= ri1) ri2++;
        int ai1 = rng_int(r, r->d12_ncnt);
        int ai2 = rng_int(r, r->d12_ncnt - 1);
        if (ai2 >= ai1) ai2++;

        int rem1 = r->d12_mem[ri1], rem2 = r->d12_mem[ri2];
        int add1 = r->d12_nmem[ai1], add2 = r->d12_nmem[ai2];

        int old_d12[MAX_M], old_d12T[MAX_M];
        memcpy(old_d12, r->delta_12, sizeof(int) * M);
        memcpy(old_d12T, r->delta_12T, sizeof(int) * M);

        r->D12[rem1] = 0; r->D12[rem2] = 0;
        r->D12[add1] = 1; r->D12[add2] = 1;
        r->D12T[mod_neg[rem1]] = 0; r->D12T[mod_neg[rem2]] = 0;
        r->D12T[mod_neg[add1]] = 1; r->D12T[mod_neg[add2]] = 1;

        compute_delta_full(r->D12, r->delta_12);
        compute_delta_full(r->D12T, r->delta_12T);

        int new_cost = compute_cost_r(r);
        int dc = new_cost - r->current_cost;

        if (dc <= 0 || rng_double(r) < exp(-(double)dc / fmax(r->temperature, 1e-10))) {
            r->current_cost = new_cost;
            d12_swap_to_nmem(r, rem1);
            d12_swap_to_nmem(r, rem2);
            d12_swap_to_mem(r, add1);
            d12_swap_to_mem(r, add2);
        } else {
            r->D12[rem1] = 1; r->D12[rem2] = 1;
            r->D12[add1] = 0; r->D12[add2] = 0;
            r->D12T[mod_neg[rem1]] = 1; r->D12T[mod_neg[rem2]] = 1;
            r->D12T[mod_neg[add1]] = 0; r->D12T[mod_neg[add2]] = 0;
            memcpy(r->delta_12, old_d12, sizeof(int) * M);
            memcpy(r->delta_12T, old_d12T, sizeof(int) * M);
        }
    }
}

/* --- Swap RNG for replica exchange --- */
static unsigned long long swap_rng_s[4];

static unsigned long long swap_rng_next(void) {
    const unsigned long long result = rotl(swap_rng_s[1] * 5, 7) * 9;
    const unsigned long long t = swap_rng_s[1] << 17;
    swap_rng_s[2] ^= swap_rng_s[0];
    swap_rng_s[3] ^= swap_rng_s[1];
    swap_rng_s[1] ^= swap_rng_s[2];
    swap_rng_s[0] ^= swap_rng_s[3];
    swap_rng_s[2] ^= t;
    swap_rng_s[3] = rotl(swap_rng_s[3], 45);
    return result;
}

static double swap_rng_double(void) {
    return (swap_rng_next() >> 11) * 0x1.0p-53;
}

/* --- Main solver with parallel tempering --- */

static int solve_pt(int d11_size, int seed, int max_iter) {
    /* Temperature ladder: geometric spacing from T_high to T_low */
    double T_high = 20.0;
    double T_low = 0.05;
    double ratio = pow(T_low / T_high, 1.0 / (NUM_REPLICAS - 1));

    /* Initialize swap RNG */
    swap_rng_s[0] = (unsigned long long)seed * 314159ULL + 271828ULL;
    swap_rng_s[1] = swap_rng_s[0] ^ 0x6C62272E07BB0142ULL;
    swap_rng_s[2] = swap_rng_s[0] ^ 0xBF58476D1CE4E5B9ULL;
    swap_rng_s[3] = swap_rng_s[0] ^ 0x94D049BB133111EBULL;
    for (int i = 0; i < 20; i++) swap_rng_next();

    /* Initialize replicas with different temperatures */
    for (int i = 0; i < NUM_REPLICAS; i++) {
        replicas[i].temperature = T_high * pow(ratio, i);
        init_replica(&replicas[i], d11_size, seed, i);
    }

    global_best_cost = 999;
    for (int i = 0; i < NUM_REPLICAS; i++) {
        if (replicas[i].current_cost < global_best_cost) {
            global_best_cost = replicas[i].current_cost;
            memcpy(best_D11, replicas[i].D11, sizeof(int) * M);
            memcpy(best_D12, replicas[i].D12, sizeof(int) * M);
        }
    }

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    int swap_accepts = 0, swap_attempts = 0;

    for (int it = 0; it < max_iter; it++) {
        /* Do one SA step on each replica */
        for (int i = 0; i < NUM_REPLICAS; i++) {
            sa_step(&replicas[i]);
        }

        /* Update global best */
        for (int i = 0; i < NUM_REPLICAS; i++) {
            if (replicas[i].current_cost < global_best_cost) {
                global_best_cost = replicas[i].current_cost;
                memcpy(best_D11, replicas[i].D11, sizeof(int) * M);
                memcpy(best_D12, replicas[i].D12, sizeof(int) * M);

                if (global_best_cost == 0) {
                    int v12_cost = verify_v1v2(&replicas[i]);
                    if (v12_cost == 0) {
                        clock_gettime(CLOCK_MONOTONIC, &ts_now);
                        double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                                       (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
                        printf("  SOLUTION at it=%d (%.1fs), replica=%d, T=%.5f\n",
                               it, elapsed, i, replicas[i].temperature);
                        return 0;
                    } else {
                        global_best_cost = v12_cost;
                    }
                }
            }
        }

        /* Attempt replica swaps every SWAP_INTERVAL iterations */
        if (it % SWAP_INTERVAL == 0 && it > 0) {
            /* Alternate even/odd pairs to ensure detailed balance */
            int parity = (it / SWAP_INTERVAL) % 2;
            for (int i = parity; i < NUM_REPLICAS - 1; i += 2) {
                swap_attempts++;
                double beta_i = 1.0 / replicas[i].temperature;
                double beta_j = 1.0 / replicas[i + 1].temperature;
                double delta_beta = beta_j - beta_i;
                double delta_cost = (double)(replicas[i].current_cost - replicas[i + 1].current_cost);

                /* Metropolis criterion for replica exchange:
                 * Accept if exp(delta_beta * delta_cost) > uniform(0,1)
                 * i.e., swap if hot replica has higher cost (favorable) */
                if (delta_beta * delta_cost >= 0 || swap_rng_double() < exp(delta_beta * delta_cost)) {
                    /* Swap temperatures */
                    double tmp_t = replicas[i].temperature;
                    replicas[i].temperature = replicas[i + 1].temperature;
                    replicas[i + 1].temperature = tmp_t;
                    swap_accepts++;
                }
            }
        }

        if (it % 1000000 == 0 && it > 0) {
            clock_gettime(CLOCK_MONOTONIC, &ts_now);
            double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                           (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;

            /* Find cost of coldest replica */
            int cold_cost = 999;
            for (int i = 0; i < NUM_REPLICAS; i++) {
                if (replicas[i].temperature < 0.06 && replicas[i].current_cost < cold_cost)
                    cold_cost = replicas[i].current_cost;
            }

            double swap_rate = swap_attempts > 0 ? (double)swap_accepts / swap_attempts : 0;
            printf("    %dM: gbest=%d cold=%d swap=%.2f (%.0f it/s)",
                   it / 1000000, global_best_cost, cold_cost, swap_rate,
                   (double)it * NUM_REPLICAS / elapsed);

            /* Print temperature of each replica compactly */
            printf(" T=[");
            for (int i = 0; i < NUM_REPLICAS; i++) {
                if (i > 0) printf(",");
                printf("%.2f", replicas[i].temperature);
            }
            printf("]");

            printf("\n");
            fflush(stdout);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double elapsed = (ts_now.tv_sec - ts_start.tv_sec) +
                   (ts_now.tv_nsec - ts_start.tv_nsec) * 1e-9;
    printf("  seed=%d: best=%d (%.1fs)\n", seed, global_best_cost, elapsed);
    fflush(stdout);

    return global_best_cost;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <n> [max_seeds] [max_iter] [fixed_d11_size]\n", argv[0]);
        return 1;
    }

    N_PARAM = atoi(argv[1]);
    int max_seeds = argc > 2 ? atoi(argv[2]) : 50;
    int max_iter = argc > 3 ? atoi(argv[3]) : 20000000;
    int fixed_d11 = argc > 4 ? atoi(argv[4]) : -1;

    M = 2 * N_PARAM - 1;
    N_TOTAL = 2 * M;
    RED_THRESH = N_PARAM - 2;
    BLUE_THRESH = N_PARAM - 1;
    D12_SIZE = N_PARAM - 1;

    printf("Parallel Tempering: n=%d, m=%d, N=%d\n", N_PARAM, M, N_TOTAL);
    printf("red_thresh=%d, blue_thresh=%d, |D12|=%d\n", RED_THRESH, BLUE_THRESH, D12_SIZE);
    printf("Replicas: %d, Seeds: %d, Iter/seed: %d\n", NUM_REPLICAS, max_seeds, max_iter);
    printf("Swap interval: %d\n", SWAP_INTERVAL);
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
            int result = solve_pt(d11_size, seed, max_iter);
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
