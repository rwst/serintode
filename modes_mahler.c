#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "solver.h"

typedef struct {
    long k;
} mahler_print_ctx;

static void print_monomial_mahler(FILE *eqs, FILE *out, long i, const void *raw_ctx)
{
    const mahler_print_ctx *ctx = (const mahler_print_ctx *) raw_ctx;
    if (i == 0L) {
        fprintf(eqs, "*f(x)");
        fprintf(out, "*f(x)");
    } else {
        long ki = 1L;
        for (long m = 0L; m < i; m++) ki *= ctx->k;
        fprintf(eqs, "*f(x^%ld)", ki);
        fprintf(out, "*f(x^%ld)", ki);
    }
}

static long compute_image_size(mpz_t *S, const char *known, long N, long *out_n_known)
{
    long unique = 0L;
    long n_known = 0L;
    for (long i = 0L; i < N; i++) {
        if (known && !known[i]) continue;
        n_known++;
        int is_new = 1;
        for (long j = 0L; j < i; j++) {
            if (known && !known[j]) continue;
            if (mpz_cmp(S[i], S[j]) == 0) { is_new = 0; break; }
        }
        if (is_new) unique++;
    }
    if (out_n_known) *out_n_known = n_known;
    return unique;
}

static const char mahler_usage[] =
    "Usage: serintode mahler <input-file> [--checks=N] [--min-order=N] [--max-coeffs=N]\n"
    "                                     [--k-min=N] [--k-max=N] [--sparse]\n";

static int mahler_run(int argc, char *argv[])
{
    long NUM_CHECKS = 6L;
    long MIN_MAHLER_ORDER = 1L;
    long MAX_COEFFS = 400L;
    long K_MIN = 2L;
    long K_MAX = 10L;
    int sparse = 0;

    static const struct option long_opts[] = {
        {"checks",     required_argument, 0, 'c'},
        {"min-order",  required_argument, 0, 'm'},
        {"max-coeffs", required_argument, 0, 'M'},
        {"k-min",      required_argument, 0, 'a'},
        {"k-max",      required_argument, 0, 'b'},
        {"sparse",     no_argument,       0, 's'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "", long_opts, NULL)) != -1) {
        switch (c) {
            case 'c': NUM_CHECKS = atol(optarg); break;
            case 'm': MIN_MAHLER_ORDER = atol(optarg); break;
            case 'M': MAX_COEFFS = atol(optarg); break;
            case 'a': K_MIN = atol(optarg); break;
            case 'b': K_MAX = atol(optarg); break;
            case 's': sparse = 1; break;
            case 'h': fputs(mahler_usage, stdout); return EXIT_SUCCESS;
            case '?': return EXIT_FAILURE;
        }
    }
    if (optind != argc - 1) {
        fputs(mahler_usage, stderr);
        return EXIT_FAILURE;
    }
    const char *finname = argv[optind];

    setvbuf(stdout, NULL, _IONBF, 0);
    time_t start, end;
    time(&start);

    mpz_t *S = NULL;
    char *known = NULL;
    long NUM_COEFFS = sparse
        ? read_sparse_series(finname, MAX_COEFFS, &S, &known)
        : read_series(finname, MAX_COEFFS, &S);

    mpz_t *M = (mpz_t *) s_malloc(MAX_COEFFS * MAX_COEFFS * sizeof(mpz_t));

    printf("File %s%s:\n", finname, sparse ? " (sparse)" : "");

    mpz_t *N = NULL;
    long nulldim = 0L;
    long mahler_order = 0L, D = 0L, k_used = 0L, COLUMNS = 0L, ROWS = 0L, numterms = 0L;
    long max_mahler_bound_for_k = 0L;
    linear_select_result picked = { -1L, 0L, 0L, 0L, 0L };
    int found = 0;

    for (long k = K_MIN; k <= K_MAX && !found; k++) {
        if (NUM_COEFFS - NUM_CHECKS <= 1L) break;
        long max_mahler_order_k = (long) floor(log((double)(NUM_COEFFS - NUM_CHECKS)) / log((double) k));
        for (long n = MIN_MAHLER_ORDER; n <= max_mahler_order_k && !found; n++) {
            long max_poly_order_n = (long) floor((double)(NUM_COEFFS - NUM_CHECKS) / (double)(n + 1L)) - 1L;
            if (max_poly_order_n < 0L) continue;

            long *ki_pow = (long *) s_malloc((n + 1L) * sizeof(long));
            ki_pow[0] = 1L;
            for (long i = 1L; i <= n; i++) ki_pow[i] = ki_pow[i - 1L] * k;

            for (long Di = 0L; Di <= max_poly_order_n; Di++) {
                numterms = n + 1L;
                COLUMNS = numterms * (Di + 1L);
                ROWS = COLUMNS + NUM_CHECKS;

                for (long col = 0L; col < COLUMNS; col++)
                    for (long row = 0L; row < ROWS; row++)
                        mpz_init(M[row*COLUMNS + col]);

                long row_max = (sparse && NUM_COEFFS > ROWS) ? NUM_COEFFS : ROWS;
                long effective_rows = 0L;
                for (long src_row = 0L; src_row < row_max && effective_rows < ROWS; src_row++) {
                    int usable = 1;
                    if (known) {
                        for (long i_term = 0L; i_term <= n && usable; i_term++) {
                            for (long j_poly = 0L; j_poly <= Di && usable; j_poly++) {
                                if (j_poly > src_row) continue;
                                long diff = src_row - j_poly;
                                if (diff % ki_pow[i_term] != 0L) continue;
                                long ai_idx = diff / ki_pow[i_term];
                                if (ai_idx >= NUM_COEFFS) continue;
                                if (!known[ai_idx]) usable = 0;
                            }
                        }
                    }
                    if (!usable) continue;
                    for (long i_term = 0L; i_term <= n; i_term++) {
                        for (long j_poly = 0L; j_poly <= Di; j_poly++) {
                            if (j_poly > src_row) continue;
                            long diff = src_row - j_poly;
                            if (diff % ki_pow[i_term] != 0L) continue;
                            long ai_idx = diff / ki_pow[i_term];
                            if (ai_idx >= NUM_COEFFS) continue;
                            long col = j_poly * numterms + i_term;
                            mpz_set(M[effective_rows*COLUMNS + col], S[ai_idx]);
                        }
                    }
                    effective_rows++;
                }

                printf("Mahler k=%ld order=%ld poly-deg=%ld checks=%ld size=%ldx%ld%s\n",
                       k, n, Di, NUM_CHECKS, effective_rows, COLUMNS,
                       (sparse && effective_rows < ROWS) ? " (sparse, fewer rows)" : "");

                if (effective_rows < COLUMNS + 1L) {
                    for (long col = 0L; col < COLUMNS; col++)
                        for (long row = 0L; row < ROWS; row++)
                            mpz_clear(M[row*COLUMNS + col]);
                    continue;
                }

                nulldim = compute_nullspace(effective_rows, COLUMNS, M, &N);

                for (long col = 0L; col < COLUMNS; col++)
                    for (long row = 0L; row < ROWS; row++)
                        mpz_clear(M[row*COLUMNS + col]);

                if (nulldim == 0L) continue;

                picked = select_best_linear(N, n, Di, nulldim, max_mahler_order_k);
                if (picked.index >= 0L) {
                    mahler_order = n;
                    D = Di;
                    k_used = k;
                    max_mahler_bound_for_k = max_mahler_order_k;
                    found = 1;
                    break;
                }
                free_nullspace(N, COLUMNS * nulldim);
                N = NULL;
            }
            free(ki_pow);
        }
    }
    (void) max_mahler_bound_for_k;

    long n_known = 0L;
    long image_size = compute_image_size(S, known, NUM_COEFFS, &n_known);
    free_series(S, NUM_COEFFS);
    free(known);
    free(M);

    if (found) {
        printf("\n***********************\n");
        printf("***Found a solution!***\n");
        printf("*Confidence level: %02ld%%*\n",
               (long) floor(100.0 - 100.0 * (mahler_order + 1L) * (picked.max_poly_order + 1L)
                            / (double)(NUM_COEFFS - mahler_order)));
        if (sparse) printf("*Image size: %ld of %ld known terms*\n", image_size, n_known);
        else        printf("*Image size: %ld of %ld terms*\n", image_size, NUM_COEFFS);
        printf("***********************\n\n");

        char fouteqsname[4096];
        snprintf(fouteqsname, sizeof(fouteqsname),
                 "%s_mahler_%ld-checks.txt", finname, NUM_CHECKS);
        FILE *fouteqs = fopen(fouteqsname, "w");
        if (fouteqs == NULL) {
            fprintf(stderr, "\nError: Could not open equations output file %s. %s\n",
                    fouteqsname, strerror(errno));
            return EXIT_FAILURE;
        }
        setvbuf(fouteqs, NULL, _IOLBF, 32);

        char header_eqs[4096];
        snprintf(header_eqs, sizeof(header_eqs), "Mahler%s := ", finname);

        mahler_print_ctx ctx = { .k = k_used };
        print_ode(fouteqs, stdout,
                  N, mahler_order + 1L, D, nulldim, picked.index,
                  header_eqs, print_monomial_mahler, &ctx);

        free_nullspace(N, COLUMNS * nulldim);
        fclose(fouteqs);
    } else {
        printf("Couldn't find any Mahler equation in k=[%ld..%ld]\n", K_MIN, K_MAX);
    }

    time(&end);
    printf("\nEllapsed time %.fs\n", difftime(end, start));
    return EXIT_SUCCESS;
}

const mode_descriptor mode_mahler = {
    .name = "mahler",
    .description = "Mahler equation search via IML (k-regular structural detector)",
    .run = mahler_run,
};
