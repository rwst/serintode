# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`serintode` searches for ODEs that annihilate an integer series (the input is a list of coefficients, one per line). It tries successively larger ODE orders / polynomial-coefficient degrees / nonlinearity depths, building a matrix whose null vectors correspond to candidate ODEs, and computing that nullspace with either IML (mod-p with verification) or FLINT (full bigint, no verification).

There are four independent main programs, each a single `.c` file with `int main()`:

- `serintode_iml.c` — linear ODE search via IML.
- `serintode_flint.c` — linear ODE search via FLINT (slower, can give spurious results without enough `NUM_CHECKS`; has no internal verification — see README).
- `serintode_iml_nonlin.c` — algebraic (nonlinear) ODE search via IML; iterates depths `1..MAX_DEPTH`. Recomputes the term-exponent table on every run.
- `serintode_iml_nonlin_lookup.c` — same as above but reads pre-generated term-exponent tables from `lookuptables/o<order>d<depth>.txt` for speed.
- `makelookup.c` — produces those lookup tables. Must be run once before `serintode_iml_nonlin_lookup`.

`old/` contains earlier prototypes — ignore unless asked.

## Build

There is no build system; compile each program directly with gcc. On this
machine IML lives under `/home/ralf/math/iml`, so the IML programs need
`-I` / `-L` overrides — the README's plain `-liml` won't find it:

```
gcc -Wall serintode_iml.c              -o serintode_iml.o              -I /home/ralf/math/iml/include -L /home/ralf/math/iml/lib64 -liml -lcblas -lgmp -lm
gcc -Wall serintode_iml_nonlin.c       -o serintode_iml_nonlin.o       -I /home/ralf/math/iml/include -L /home/ralf/math/iml/lib64 -liml -lcblas -lgmp -lm
gcc -Wall serintode_iml_nonlin_lookup.c -o serintode_iml_nonlin_lookup.o -I /home/ralf/math/iml/include -L /home/ralf/math/iml/lib64 -liml -lcblas -lgmp -lm
gcc -Wall makelookup.c                 -o makelookup                   -lgmp -lm
```

The `.o` extension is misleading — these are full executables. Library install instructions (IML, GMP, ATLAS) are in `README.txt`.

## How to run

All in-scope programs take the input file path as `argv[1]`:

```
./serintode_iml.o              <input-file>
./serintode_iml_nonlin.o       <input-file>
./serintode_iml_nonlin_lookup.o <input-file>
./makelookup                   <MAX_ODE_ORDER> <NUM_COEFFS>     # writes lookuptables/o<n>d<p>.txt
```

The other tunables — `NUM_CHECKS`, `MIN_ODE_ORDER`, `MAX_COEFFS`, `MAX_LINE_LENGTH`, plus `MIN_DEPTH`/`MAX_DEPTH` on the nonlin variants — are still `const` locals at the top of `main()` and must be edited in the source, then recompiled. The relevant block in each file:

```c
long const NUM_CHECKS = ...;
long const MIN_ODE_ORDER = ...;
long const MAX_COEFFS = ...;
long const MAX_LINE_LENGTH = 100000L;
// nonlin variants additionally have MIN_DEPTH, MAX_DEPTH
```

`tests/central_binomials.txt` is a known-good regression input (OEIS A000984; expected ODE `(1-4x)y' - 2y = 0`). Input files are plain text, one base-10 integer per line; leading zero coefficients are stripped automatically.

Output: when a solution is found, the program writes `<finname>_solution_<NUM_CHECKS>-checks.txt` (linear) or `<finname>_nonlinsol_<NUM_CHECKS>-checks.txt` (nonlin) in Maple syntax, alongside printing it to stdout.

There is no test suite.

## Architecture notes worth knowing before editing

**The main loop structure is the same in all four programs**: outer loop over ODE order `n`, inner loop builds a matrix `M` whose rows encode the ODE-applied-to-series constraints, computes its right nullspace (`kernelMP`/`nullspaceMP` from IML or `fmpz_mat_nullspace` from FLINT), and if the nullspace is nontrivial picks a "best" null vector (by minimising max polynomial degree, then number of nonzero terms) as the candidate ODE.

**The IML constraint** (`README.txt` lines 80–83): IML reduces mod a random 15–19-bit prime and uses CBLAS doubles internally, which caps `ceil(n/2)*(p-1)^2 + (p-1) ≤ 2^53 - 1`. This is why `MAX_COEFFS` matters — too large and IML's row-echelon step loses precision.

**FLINT version is known-imperfect**: README and inline comments warn it can return spurious results unless `NUM_CHECKS` is large enough, because there's no post-construction annihilation check. IML's output is self-verified.

**Nonlin term enumeration**: both nonlin programs use the recursive `combs()` function (credited to Richard Brent) to enumerate compositions — the multi-indices for terms like `y^a * (Dy)^b * (D²y)^c`. The `_lookup` variant precomputes these into files; the non-lookup variant regenerates them every iteration.

**MAX_POLY_ORDER formula**: `floor((NUM_COEFFS - NUM_CHECKS - ODE_ORDER) / (ODE_ORDER + 1)) - 1`. **MAX_ODE_ORDER**: `floor((NUM_COEFFS - NUM_CHECKS) / 2) - 1`. These are derived in the README.

## Known issues / gotchas

See `TODO.md` for a list of robustness items found in an input-parsing audit (unchecked GMP/FLINT parser returns, latent `sprintf`-into-64-byte-buffer overflows on long `finname`, bitwise-`&` typo in `makelookup.c`, etc.). None are reachable from the data file.

When changing parameters at the top of `main()`, remember:
- `MAX_LINE_LENGTH` is allocated as a stack array (`char input_string[MAX_LINE_LENGTH+1]`) — keep it well under the stack limit.
- `MAX_COEFFS` drives an `MAX_COEFFS*MAX_COEFFS` `mpz_t` matrix allocation, which is large fast.
- The four programs do not share code; a fix often needs to be applied in all four (and `old/` if relevant).
