#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "solver.h"

static const char kkernel_usage[] =
    "Usage: serintode kkernel <input-file> [--k=N] [--max-depth=N] [--max-coeffs=N] [--sparse]\n";

static int kkernel_run(int argc, char *argv[])
{
    long K = 2L;
    long MAX_DEPTH = 6L;
    long MAX_COEFFS = 1000L;
    int sparse = 0;

    static const struct option long_opts[] = {
        {"k",          required_argument, 0, 'k'},
        {"max-depth",  required_argument, 0, 'd'},
        {"max-coeffs", required_argument, 0, 'M'},
        {"sparse",     no_argument,       0, 's'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "", long_opts, NULL)) != -1) {
        switch (c) {
            case 'k': K = atol(optarg); break;
            case 'd': MAX_DEPTH = atol(optarg); break;
            case 'M': MAX_COEFFS = atol(optarg); break;
            case 's': sparse = 1; break;
            case 'h': fputs(kkernel_usage, stdout); return EXIT_SUCCESS;
            case '?': return EXIT_FAILURE;
        }
    }
    if (optind != argc - 1) {
        fputs(kkernel_usage, stderr);
        return EXIT_FAILURE;
    }
    if (K < 2L) {
        fprintf(stderr, "Error: --k must be >= 2.\n");
        return EXIT_FAILURE;
    }
    const char *finname = argv[optind];

    setvbuf(stdout, NULL, _IONBF, 0);

    mpz_t *S = NULL;
    char *known = NULL;
    long N = sparse
        ? read_sparse_series(finname, MAX_COEFFS, &S, &known)
        : read_series(finname, MAX_COEFFS, &S);
    printf("File %s%s, k=%ld, depths 0..%ld, N=%ld\n",
           finname, sparse ? " (sparse)" : "", K, MAX_DEPTH, N);
    printf("(rank of Q-span of {a(k^i*n + j) : 0<=i<=depth, 0<=j<k^i})\n");
    printf("k-regular  <=>  rank stabilises as depth grows.\n\n");

    for (long depth = 0L; depth <= MAX_DEPTH; depth++) {
        long total_rows = 0L;
        long ki = 1L;
        for (long i = 0L; i <= depth; i++) {
            total_rows += ki;
            ki *= K;
        }

        long L = N;
        long power = 1L;
        for (long i = 0L; i <= depth; i++) {
            for (long j = 0L; j < power; j++) {
                long row_len = (j < N) ? (N - j + power - 1L) / power : 0L;
                if (row_len < L) L = row_len;
            }
            power *= K;
        }
        if (L <= 0L) {
            printf("  depth=%ld: N=%ld too small (k^%ld=%ld)\n",
                   depth, N, depth, power / K);
            continue;
        }

        long L_eff = L;
        char *col_usable = NULL;
        if (known) {
            col_usable = (char *) s_malloc(L * sizeof(char));
            for (long n = 0L; n < L; n++) col_usable[n] = 1;
            long pwr = 1L;
            for (long i = 0L; i <= depth; i++) {
                for (long j = 0L; j < pwr; j++) {
                    for (long n = 0L; n < L; n++) {
                        long src = pwr * n + j;
                        if (src >= N || !known[src]) col_usable[n] = 0;
                    }
                }
                pwr *= K;
            }
            L_eff = 0L;
            for (long n = 0L; n < L; n++) if (col_usable[n]) L_eff++;
        }
        if (L_eff <= 0L) {
            printf("  depth=%ld: rows=%ld, common-length=%ld, no usable columns\n",
                   depth, total_rows, L);
            free(col_usable);
            continue;
        }

        long rows = total_rows;
        mpz_t *M = (mpz_t *) s_malloc(rows * L_eff * sizeof(mpz_t));
        for (long i = 0L; i < rows*L_eff; i++) mpz_init(M[i]);

        long row_idx = 0L;
        power = 1L;
        for (long i = 0L; i <= depth; i++) {
            for (long j = 0L; j < power; j++) {
                long col_dst = 0L;
                for (long n = 0L; n < L; n++) {
                    if (col_usable && !col_usable[n]) continue;
                    long src = power * n + j;
                    if (src < N) mpz_set(M[row_idx*L_eff + col_dst], S[src]);
                    col_dst++;
                }
                row_idx++;
            }
            power *= K;
        }

        long rank = compute_rank(rows, L_eff, M);

        for (long i = 0L; i < rows*L_eff; i++) mpz_clear(M[i]);
        free(M);
        free(col_usable);

        if (sparse && L_eff < L) {
            printf("  depth=%ld: rows=%ld, common-length=%ld (of %ld), rank=%ld\n",
                   depth, rows, L_eff, L, rank);
        } else {
            printf("  depth=%ld: rows=%ld, common-length=%ld, rank=%ld\n",
                   depth, rows, L, rank);
        }
    }

    free_series(S, N);
    free(known);
    return EXIT_SUCCESS;
}

const mode_descriptor mode_kkernel = {
    .name = "kkernel",
    .description = "k-kernel rank diagnostic (k-regular if rank stabilises)",
    .run = kkernel_run,
};
