#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gmp.h"

#include "serintode.h"
#include "io.h"
#include "solver.h"

static void print_monomial_linear(FILE *eqs, FILE *out, long i, const void *ctx)
{
    (void) ctx;
    if (i == 0L) {
        fprintf(eqs, "*y(x)");
        fprintf(out, "*y");
    } else if (i == 1L) {
        fprintf(eqs, "*diff(y(x),x)");
        fprintf(out, "*Dx");
    } else {
        fprintf(eqs, "*diff(y(x),x$%ld)", i);
        fprintf(out, "*Dx^%ld", i);
    }
}

static int linear_run(int argc, char *argv[])
{
    long const NUM_CHECKS = 6L;
    long const MIN_ODE_ORDER = 1L;
    long const MAX_COEFFS = 400L;

    setvbuf(stdout, NULL, _IONBF, 0);
    time_t start, end;
    time(&start);

    if (argc != 2) {
        fprintf(stderr, "Usage: serintode linear <input-file>\n");
        return EXIT_FAILURE;
    }
    const char *finname = argv[1];

    mpz_t *S = NULL;
    long NUM_COEFFS = read_series(finname, MAX_COEFFS, &S);

    mpz_t *M = (mpz_t *) s_malloc(MAX_COEFFS * MAX_COEFFS * sizeof(mpz_t));

    printf("File %s:\n", finname);
    long MAX_ODE_ORDER = (long) floor((double)(NUM_COEFFS - NUM_CHECKS) / 2L) - 1L;

    mpz_t *N = NULL;
    long nulldim = 0L;
    long ODE_ORDER = 0L;
    long MAX_POLY_ORDER = 0L;
    long COLUMNS = 0L, ROWS = 0L;
    linear_select_result picked = { -1L, 0L, 0L, 0L, 0L };

    for (long n = MIN_ODE_ORDER; n < MAX_ODE_ORDER + 1L; n++) {
        ODE_ORDER = n;
        MAX_POLY_ORDER = (long) floor((double)(NUM_COEFFS - NUM_CHECKS - ODE_ORDER) / (ODE_ORDER + 1L)) - 1L;
        if (MAX_POLY_ORDER == 0L) break;

        COLUMNS = (ODE_ORDER + 1L) * (MAX_POLY_ORDER + 1L);
        ROWS = COLUMNS + NUM_CHECKS;
        printf("Checking for ODE order %ld with polynom coeffs of order %ld and %ld checks\n",
               ODE_ORDER, MAX_POLY_ORDER, NUM_CHECKS);

        for (long i = 0L; i < COLUMNS; i++)
            for (long j = 0L; j < ROWS; j++)
                mpz_init(M[j*COLUMNS + i]);

        mpz_t temp, temp2, coeff;
        mpz_inits(temp, temp2, coeff, NULL);
        for (long i = 0L; i < MAX_POLY_ORDER + 1L; i++) {
            for (long j = 0L; j < ODE_ORDER + 1L; j++) {
                for (long k = i; k < ROWS; k++) {
                    mpz_fac_ui(temp, j + k - i);
                    mpz_fac_ui(temp2, k - i);
                    mpz_divexact(coeff, temp, temp2);
                    mpz_mul(temp, coeff, S[j + k - i]);
                    mpz_set(M[k*COLUMNS + i*(ODE_ORDER + 1L) + j], temp);
                }
            }
        }
        mpz_clears(temp, temp2, coeff, NULL);

        nulldim = compute_nullspace(ROWS, COLUMNS, M, &N);

        for (long i = 0L; i < COLUMNS; i++)
            for (long j = 0L; j < ROWS; j++)
                mpz_clear(M[j*COLUMNS + i]);

        if (nulldim == 0L) continue;

        picked = select_best_linear(N, ODE_ORDER, MAX_POLY_ORDER, nulldim, MAX_ODE_ORDER);
        if (picked.index >= 0L) break;

        free_nullspace(N, COLUMNS * nulldim);
        N = NULL;
    }

    free_series(S, NUM_COEFFS);
    free(M);

    if (picked.index >= 0L) {
        printf("\n***********************\n");
        printf("***Found a solution!***\n");
        printf("*Confidence level: %02ld%%*\n",
               (long) floor(100.0 - 100.0 * (ODE_ORDER + 1L) * (picked.max_poly_order + 1L)
                            / (double)(NUM_COEFFS - ODE_ORDER)));
        printf("***********************\n\n");

        char fouteqsname[4096];
        snprintf(fouteqsname, sizeof(fouteqsname),
                 "%s_linear_%ld-checks.txt", finname, NUM_CHECKS);
        FILE *fouteqs = fopen(fouteqsname, "w");
        if (fouteqs == NULL) {
            fprintf(stderr, "\nError: Could not open equations output file %s. %s\n",
                    fouteqsname, strerror(errno));
            return EXIT_FAILURE;
        }
        setvbuf(fouteqs, NULL, _IOLBF, 32);

        char header_eqs[4096];
        snprintf(header_eqs, sizeof(header_eqs), "ODE%s := ", finname);

        print_ode(fouteqs, stdout,
                  N, ODE_ORDER + 1L, MAX_POLY_ORDER, nulldim, picked.index,
                  header_eqs, print_monomial_linear, NULL);

        free_nullspace(N, COLUMNS * nulldim);
        fclose(fouteqs);
    } else {
        if (ODE_ORDER > 0L) printf("Couldn't find any solution\n");
        else printf("Series too short %ld checks\n", NUM_CHECKS);
    }

    time(&end);
    printf("\nEllapsed time %.fs\n", difftime(end, start));
    return EXIT_SUCCESS;
}

const mode_descriptor mode_linear = {
    .name = "linear",
    .description = "Linear ODE search via IML",
    .run = linear_run,
};
