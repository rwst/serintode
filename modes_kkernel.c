#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "solver.h"

static const char kkernel_usage[] =
    "Usage: serintode kkernel <input-file> [--k=N] [--max-depth=N] [--max-coeffs=N]\n";

static int kkernel_run(int argc, char *argv[])
{
    long K = 2L;
    long MAX_DEPTH = 6L;
    long MAX_COEFFS = 1000L;

    static const struct option long_opts[] = {
        {"k",          required_argument, 0, 'k'},
        {"max-depth",  required_argument, 0, 'd'},
        {"max-coeffs", required_argument, 0, 'M'},
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
    long N = read_series(finname, MAX_COEFFS, &S);
    printf("File %s, k=%ld, depths 0..%ld, N=%ld\n", finname, K, MAX_DEPTH, N);
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

        long rows = total_rows;
        mpz_t *M = (mpz_t *) s_malloc(rows * L * sizeof(mpz_t));
        for (long i = 0L; i < rows*L; i++) mpz_init(M[i]);

        long row_idx = 0L;
        power = 1L;
        for (long i = 0L; i <= depth; i++) {
            for (long j = 0L; j < power; j++) {
                for (long n = 0L; n < L; n++) {
                    long src = power * n + j;
                    if (src < N) mpz_set(M[row_idx*L + n], S[src]);
                }
                row_idx++;
            }
            power *= K;
        }

        long rank = compute_rank(rows, L, M);

        for (long i = 0L; i < rows*L; i++) mpz_clear(M[i]);
        free(M);

        printf("  depth=%ld: rows=%ld, common-length=%ld, rank=%ld\n",
               depth, rows, L, rank);
    }

    free_series(S, N);
    return EXIT_SUCCESS;
}

const mode_descriptor mode_kkernel = {
    .name = "kkernel",
    .description = "k-kernel rank diagnostic (k-regular if rank stabilises)",
    .run = kkernel_run,
};
