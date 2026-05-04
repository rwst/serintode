#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "solver.h"
#include "modes_nonlin.h"

/* Thanks to Richard Brent (Newcastle University) for this composition
 * enumerator. Fills r[0..C-1] with the C = (k+p)!/(k!p!) tuples of
 * length k whose entries sum to <= p (in lex-decreasing order). */
void combs(long **r, long *s, long k, long p, long q, long *arrindex)
{
    if (k <= 0) {
        for (long m = 0; m < q; m++) r[*arrindex][m] = s[m];
        (*arrindex)++;
        return;
    }
    for (long j = p; j >= 0; j--) {
        s[q] = j;
        combs(r, s, k - 1, p - j, q + 1, arrindex);
    }
}

typedef struct {
    long ode_order;
    long **orderexp;
} nonlin_print_ctx;

static void print_monomial_nonlin(FILE *eqs, FILE *out, long i, const void *raw_ctx)
{
    const nonlin_print_ctx *ctx = (const nonlin_print_ctx *) raw_ctx;
    for (long j = 0L; j < ctx->ode_order + 1L; j++) {
        long e = ctx->orderexp[i][j];
        if (e == 0L) continue;
        if (j == 0L) {
            if (e != 1L) {
                fprintf(eqs, "*(y(x))^%ld", e);
                fprintf(out, "*y^%ld", e);
            } else {
                fprintf(eqs, "*y(x)");
                fprintf(out, "*y");
            }
        } else if (j == 1L) {
            if (e != 1L) {
                fprintf(eqs, "*(diff(y(x),x))^%ld", e);
                fprintf(out, "*(Dx)^%ld", e);
            } else {
                fprintf(eqs, "*diff(y(x),x)");
                fprintf(out, "*(Dx)");
            }
        } else {
            if (e != 1L) {
                fprintf(eqs, "*(diff(y(x),x$%ld))^%ld", j, e);
                fprintf(out, "*(Dx^%ld)^%ld", j, e);
            } else {
                fprintf(eqs, "*diff(y(x),x$%ld)", j);
                fprintf(out, "*(Dx^%ld)", j);
            }
        }
    }
}

static int load_orderexp_from_lookup(const char *lookup_dir,
    long ode_order, long depth, long numterms, long **orderexp)
{
    char fname[4096];
    snprintf(fname, sizeof(fname), "%s/o%ldd%ld.txt", lookup_dir, ode_order, depth);
    FILE *flook = fopen(fname, "r");
    if (flook == NULL) {
        printf("No lookup table for order %ld depth %ld. %s\n",
               ode_order, depth, strerror(errno));
        printf("Skipping ...\n");
        return -1;
    }
    for (long i = 0L; i < numterms; i++) {
        for (long j = 0L; j <= ode_order; j++) {
            if (fscanf(flook, "%ld ", &orderexp[i][j]) != 1) {
                fprintf(stderr, "Error reading lookup table %s at i=%ld j=%ld.\n",
                        fname, i, j);
                fclose(flook);
                return -1;
            }
        }
    }
    fclose(flook);
    return 0;
}

static const char nonlin_usage[] =
    "Usage: serintode nonlin <input-file> [--checks=N] [--min-order=N] [--max-coeffs=N]\n"
    "                                     [--min-depth=N] [--max-depth=N] [--lookup-dir=PATH]\n";

static int nonlin_run(int argc, char *argv[])
{
    long NUM_CHECKS = 0L;
    long MIN_ODE_ORDER = 1L;
    long MIN_DEPTH = 1L;
    long MAX_DEPTH = 10L;
    long MAX_COEFFS = 100L;
    const char *lookup_dir = NULL;

    static const struct option long_opts[] = {
        {"checks",     required_argument, 0, 'c'},
        {"min-order",  required_argument, 0, 'm'},
        {"max-coeffs", required_argument, 0, 'M'},
        {"min-depth",  required_argument, 0, 'd'},
        {"max-depth",  required_argument, 0, 'D'},
        {"lookup-dir", required_argument, 0, 'L'},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    optind = 1;
    int c;
    while ((c = getopt_long(argc, argv, "", long_opts, NULL)) != -1) {
        switch (c) {
            case 'c': NUM_CHECKS = atol(optarg); break;
            case 'm': MIN_ODE_ORDER = atol(optarg); break;
            case 'M': MAX_COEFFS = atol(optarg); break;
            case 'd': MIN_DEPTH = atol(optarg); break;
            case 'D': MAX_DEPTH = atol(optarg); break;
            case 'L': lookup_dir = optarg; break;
            case 'h': fputs(nonlin_usage, stdout); return EXIT_SUCCESS;
            case '?': return EXIT_FAILURE;
        }
    }
    if (optind != argc - 1) {
        fputs(nonlin_usage, stderr);
        return EXIT_FAILURE;
    }
    const char *finname = argv[optind];

    setvbuf(stdout, NULL, _IONBF, 0);
    time_t start, end;
    time(&start);

    mpz_t *I = NULL;
    long NUM_COEFFS = read_series(finname, MAX_COEFFS, &I);
    printf("\nChecking %s:\n", finname);

    mpz_t *M = (mpz_t *) s_malloc(MAX_COEFFS * MAX_COEFFS * sizeof(mpz_t));
    mpz_t *temparray = (mpz_t *) s_malloc(MAX_COEFFS * sizeof(mpz_t));
    for (long i = 0L; i < MAX_COEFFS; i++) mpz_init(temparray[i]);

    long MAX_ODE_ORDER = NUM_COEFFS - NUM_CHECKS - 1L;
    long maxnumterms = NUM_COEFFS - NUM_CHECKS;

    /* D[i][k] = i-th derivative coefficient at index k. */
    mpz_t **D = (mpz_t **) s_malloc((MAX_ODE_ORDER + 1L) * sizeof(mpz_t *));
    for (long i = 0L; i < MAX_ODE_ORDER + 1L; i++) {
        D[i] = (mpz_t *) s_calloc(NUM_COEFFS, sizeof(mpz_t));
        for (long j = 0L; j < NUM_COEFFS; j++) mpz_init(D[i][j]);
    }
    for (long i = 0L; i < NUM_COEFFS; i++) mpz_set(D[0L][i], I[i]);
    {
        mpz_t temp;
        mpz_init(temp);
        for (long i = 1L; i < MAX_ODE_ORDER + 1L; i++) {
            for (long j = 0L; j < NUM_COEFFS - i; j++) {
                mpz_mul_ui(temp, D[i - 1L][j + 1L], (unsigned long)(j + 1L));
                mpz_set(D[i][j], temp);
            }
        }
        mpz_clear(temp);
    }

    mpz_t *N = NULL;
    long nulldim = 0L;
    long ODE_ORDER = 0L;
    long MAX_POLY_ORDER = 0L;
    long COLUMNS = 0L, ROWS = 0L;
    long numterms = 0L;
    long depth_used = 0L;
    long **orderexp_kept = NULL;
    long orderexp_kept_rows = 0L;
    long *s_buf = NULL;
    nonlin_select_result picked = { -1L, 0L, 0L, 0L, 0L, 0L };

    for (long n = MIN_ODE_ORDER; n < MAX_ODE_ORDER + 1L; n++) {
        ODE_ORDER = n;

        long MAX_DEPTH_POSS = 1L;
        long ordermaxnumterms = 0L;
        {
            mpz_t temp, temp2;
            mpz_inits(temp, temp2, NULL);
            for (MAX_DEPTH_POSS = 1L; MAX_DEPTH_POSS <= maxnumterms; MAX_DEPTH_POSS++) {
                mpz_fac_ui(temp, ODE_ORDER + 1L + MAX_DEPTH_POSS);
                mpz_fac_ui(temp2, ODE_ORDER + 1L);
                mpz_divexact(temp, temp, temp2);
                mpz_fac_ui(temp2, MAX_DEPTH_POSS);
                mpz_divexact(temp, temp, temp2);
                if ((long) mpz_get_ui(temp) > maxnumterms) {
                    MAX_DEPTH_POSS--;
                    break;
                }
                ordermaxnumterms = (long) mpz_get_ui(temp);
            }
            mpz_clears(temp, temp2, NULL);
        }
        if (MAX_DEPTH_POSS > MAX_DEPTH) MAX_DEPTH_POSS = MAX_DEPTH;

        mpz_t **S = (mpz_t **) s_malloc(ordermaxnumterms * sizeof(mpz_t *));
        for (long i = 0L; i < ordermaxnumterms; i++) {
            S[i] = (mpz_t *) s_calloc(NUM_COEFFS, sizeof(mpz_t));
            for (long j = 0L; j < NUM_COEFFS; j++) mpz_init(S[i][j]);
        }

        int found = 0;
        for (long p = MIN_DEPTH; p < MAX_DEPTH_POSS + 1L; p++) {
            mpz_t temp, temp2;
            mpz_inits(temp, temp2, NULL);
            mpz_fac_ui(temp, ODE_ORDER + 1L + p);
            mpz_fac_ui(temp2, ODE_ORDER + 1L);
            mpz_divexact(temp, temp, temp2);
            mpz_fac_ui(temp2, p);
            mpz_divexact(temp, temp, temp2);
            numterms = (long) mpz_get_ui(temp) - 1L;
            mpz_clears(temp, temp2, NULL);

            MAX_POLY_ORDER = (long) floor((double)(NUM_COEFFS - NUM_CHECKS - ODE_ORDER) / numterms) - 1L;
            if (MAX_POLY_ORDER < 1L) break;

            long *s = (long *) s_calloc(ODE_ORDER + 1L, sizeof(long));
            long **orderexp = (long **) s_malloc((numterms + 1L) * sizeof(long *));
            for (long i = 0L; i < numterms + 1L; i++) {
                orderexp[i] = (long *) s_calloc(ODE_ORDER + 1L, sizeof(long));
            }

            if (lookup_dir != NULL) {
                if (load_orderexp_from_lookup(lookup_dir, ODE_ORDER, p, numterms, orderexp) != 0) {
                    for (long i = 0L; i < numterms + 1L; i++) free(orderexp[i]);
                    free(orderexp);
                    free(s);
                    break;
                }
            } else {
                long arrindex = 0L;
                combs(orderexp, s, ODE_ORDER + 1L, p, 0L, &arrindex);
            }

            /* Build S[i][k] = coeff of x^k in the i-th nonlinear monomial. */
            for (long i = 0L; i < numterms; i++) {
                int first = 0;
                for (long j = 0L; j < ODE_ORDER + 1L; j++) {
                    if (orderexp[i][j] == 0L) continue;
                    if (first == 0) {
                        for (long k = 0L; k < NUM_COEFFS - j; k++) mpz_set(S[i][k], D[j][k]);
                        first++;
                        for (long m = 1L; m < orderexp[i][j]; m++) {
                            for (long k = 0L; k < NUM_COEFFS - j; k++) {
                                mpz_t acc, prod;
                                mpz_inits(acc, prod, NULL);
                                for (long l = 0L; l < k + 1L; l++) {
                                    mpz_mul(prod, S[i][l], D[j][k - l]);
                                    mpz_add(acc, acc, prod);
                                }
                                mpz_set(temparray[k], acc);
                                mpz_clears(acc, prod, NULL);
                            }
                            for (long k = 0L; k < NUM_COEFFS - j; k++) mpz_set(S[i][k], temparray[k]);
                        }
                    } else {
                        for (long m = 0L; m < orderexp[i][j]; m++) {
                            for (long k = 0L; k < NUM_COEFFS; k++) {
                                mpz_t acc, prod;
                                mpz_inits(acc, prod, NULL);
                                for (long l = 0L; l < k + 1L; l++) {
                                    mpz_mul(prod, S[i][l], D[j][k - l]);
                                    mpz_add(acc, acc, prod);
                                }
                                mpz_set(temparray[k], acc);
                                mpz_clears(acc, prod, NULL);
                            }
                            for (long k = 0L; k < NUM_COEFFS - j; k++) mpz_set(S[i][k], temparray[k]);
                        }
                    }
                }
            }

            COLUMNS = numterms * (MAX_POLY_ORDER + 1L);
            ROWS = COLUMNS + NUM_CHECKS;
            printf("ODE of order %ld with %ld nonlin terms depth %ld, polynom coeffs of order %ld and %ld checks, size %ld x %ld\n",
                   ODE_ORDER, numterms, p, MAX_POLY_ORDER, NUM_CHECKS, ROWS, COLUMNS);

            for (long i = 0L; i < COLUMNS; i++)
                for (long j = 0L; j < ROWS; j++)
                    mpz_init(M[j*COLUMNS + i]);

            for (long i = 0L; i < MAX_POLY_ORDER + 1L; i++) {
                for (long j = i; j < ROWS; j++) {
                    for (long k = 0L; k < numterms; k++) {
                        mpz_set(M[j*COLUMNS + i*numterms + k], S[k][j - i]);
                    }
                }
            }

            nulldim = compute_nullspace(ROWS, COLUMNS, M, &N);

            for (long i = 0L; i < COLUMNS; i++)
                for (long j = 0L; j < ROWS; j++)
                    mpz_clear(M[j*COLUMNS + i]);

            if (nulldim > 0L) {
                picked = select_best_nonlin(N, orderexp, ODE_ORDER, numterms, MAX_POLY_ORDER, MAX_DEPTH, nulldim);
                if (picked.index >= 0L) {
                    found = 1;
                    depth_used = p;
                    orderexp_kept = orderexp;
                    orderexp_kept_rows = numterms + 1L;
                    s_buf = s;
                    break;
                }
                free_nullspace(N, COLUMNS * nulldim);
                N = NULL;
            }

            for (long i = 0L; i < numterms + 1L; i++) free(orderexp[i]);
            free(orderexp);
            free(s);
        }

        for (long i = 0L; i < ordermaxnumterms; i++) {
            for (long j = 0L; j < NUM_COEFFS; j++) mpz_clear(S[i][j]);
            free(S[i]);
        }
        free(S);

        if (found) break;
    }

    for (long i = 0L; i < MAX_ODE_ORDER + 1L; i++) {
        for (long j = 0L; j < NUM_COEFFS; j++) mpz_clear(D[i][j]);
        free(D[i]);
    }
    free(D);
    free_series(I, NUM_COEFFS);
    for (long i = 0L; i < MAX_COEFFS; i++) mpz_clear(temparray[i]);
    free(temparray);
    free(M);

    if (picked.index >= 0L) {
        printf("null dimension = %ld\n", nulldim);
        printf("****Found a solution!****\n");

        char fouteqsname[4096];
        snprintf(fouteqsname, sizeof(fouteqsname),
                 "%s_nonlin_%ld-checks.txt", finname, NUM_CHECKS);
        FILE *fouteqs = fopen(fouteqsname, "w");
        if (fouteqs == NULL) {
            fprintf(stderr, "\nError: Could not open equations output file %s. %s\n",
                    fouteqsname, strerror(errno));
            return EXIT_FAILURE;
        }
        setvbuf(fouteqs, NULL, _IOLBF, 32);

        char header_eqs[4096];
        snprintf(header_eqs, sizeof(header_eqs), "#ODE candidate for %s\n", finname);

        nonlin_print_ctx ctx = { .ode_order = ODE_ORDER, .orderexp = orderexp_kept };
        print_ode(fouteqs, stdout,
                  N, numterms, MAX_POLY_ORDER, nulldim, picked.index,
                  header_eqs, print_monomial_nonlin, &ctx);

        printf("Confidence level: %02ld%%\n",
               (long) floor(100.0 - 100.0 * picked.terms_used / (double)(NUM_COEFFS - ODE_ORDER)));

        free_nullspace(N, COLUMNS * nulldim);
        fclose(fouteqs);

        if (orderexp_kept != NULL) {
            for (long i = 0L; i < orderexp_kept_rows; i++) free(orderexp_kept[i]);
            free(orderexp_kept);
            free(s_buf);
        }
        (void) depth_used;
    } else {
        if (ODE_ORDER > 0L) printf("Couldn't find any solution\n");
        else printf("Series too short %ld checks\n", NUM_CHECKS);
    }

    time(&end);
    printf("\nEllapsed time %.fs\n", difftime(end, start));
    return EXIT_SUCCESS;
}

const mode_descriptor mode_nonlin = {
    .name = "nonlin",
    .description = "Algebraic (nonlinear) ODE search via IML",
    .run = nonlin_run,
};
