# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`serintode` searches for ODEs that annihilate an integer series (the input is a list of coefficients, one per line). It tries successively larger ODE orders / polynomial-coefficient degrees / nonlinearity depths, building a matrix whose null vectors correspond to candidate ODEs, and computing that nullspace with either IML (mod-p with verification) or FLINT (full bigint, no verification).

A staged refactor (see `plan.md`) is in progress, collapsing the IML programs into a single dispatcher binary `serintode` with subcommands. After Stage 4:

- `serintode linear <input>` — linear ODE search via IML. In `modes_linear.c`.
- `serintode nonlin <input> [--lookup-dir=PATH]` — algebraic (nonlinear) ODE search via IML. In `modes_nonlin.c`. With `--lookup-dir` set, reads pre-generated term-exponent tables from `<PATH>/o<n>d<p>.txt`; without it, enumerates compositions in-process via `combs()`.
- `serintode makelookup <max-order> <num-coeffs>` — generates the lookup tables. In `modes_makelookup.c`, sharing `combs()` from `modes_nonlin.h`.
- Shared code: `io.{c,h}` (parser, pretty-printer, `s_malloc`), `solver.{c,h}` (kernel wrapper, null-vector selectors).
- `serintode_flint.c` — out of scope (FLINT path).
- Orphaned, pending removal: `serintode_iml.c`, `serintode_iml_nonlin.c`, `serintode_iml_nonlin_lookup.c`, `makelookup.c`. Stage 5 deletes them.

`old/` contains earlier prototypes — ignore unless asked.

## Build

```
make            # builds all four programs
make test       # runs the regression suite
make clean
```

Override `IML_DIR=/path/to/iml` if IML lives somewhere other than `/home/ralf/math/iml`. The `.o` extension on the binaries is misleading — these are full executables. Library install instructions (IML, GMP, ATLAS) are in `README.txt`.

## How to run

```
./serintode linear     <input-file>
./serintode nonlin     <input-file> [--lookup-dir=PATH]
./serintode makelookup <MAX_ODE_ORDER> <NUM_COEFFS>     # writes lookuptables/o<n>d<p>.txt
```

Tunables (`NUM_CHECKS`, `MIN_ODE_ORDER`, `MAX_COEFFS`, plus `MIN_DEPTH`/`MAX_DEPTH` on nonlin) are still `const` locals at the top of each mode's run function and must be edited in the source, then recompiled. `MAX_LINE_LENGTH` is a `#define` in `io.h`. Stage 5 escalates the per-run knobs to `getopt_long`.

`tests/central_binomials.txt` is a known-good regression input (OEIS A000984; expected ODE `(1-4x)y' - 2y = 0`). Input files are plain text, one base-10 integer per line; leading zero coefficients are stripped automatically.

Output filename: `<finname>_<mode>_<NUM_CHECKS>-checks.txt` in Maple syntax, also echoed to stdout.

`make test` runs `./serintode linear` and `./serintode nonlin` against `tests/central_binomials.txt`, diffing against `tests/expected/`.

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
