# Refactor plan: collapse the four serintode_* programs into one

Resolves the **Maintainability** item in `TODO.md`. Lands before the new
`modes_mahler.c` so the Mahler port doesn't become the fourth offender.

## Goals

- One binary with mode dispatch instead of three near-duplicate IML
  programs plus `makelookup`.
- Shared code for the input parser, the polynomial-coefficient
  pretty-printer, the IML kernel call, the spurious-zero check, and the
  null-vector selection harness.
- A regression net (`make test`) so the rewrite can be verified
  end-to-end, not just by eyeballing.
- A clean place for `modes_mahler.c` to land in a follow-up PR.

In scope: `serintode_iml.c`, `serintode_iml_nonlin.c`,
`serintode_iml_nonlin_lookup.c`, `makelookup.c`. Out of scope:
`serintode_flint.c`, `old/`.

## Architecture decisions

### Shape: one binary, subcommand dispatch

```
serintode linear      <input>    [--checks=N] [--min-order=N] [--max-coeffs=N]
serintode nonlin      <input>    [--checks=N] [--min-order=N] [--max-coeffs=N]
                                 [--min-depth=N] [--max-depth=N] [--lookup-dir=PATH]
serintode mahler      <input>    [--checks=N] [--min-order=N] [--max-coeffs=N]
                                 [--k-min=N] [--k-max=N]
serintode makelookup  <max-order> <num-coeffs>
```

The two existing `nonlin` programs collapse into one `nonlin` mode.
`--lookup-dir=PATH` reads `<PATH>/o<n>d<p>.txt`; absent, `combs()`
enumerates compositions in-process.

### File layout (6 files)

```
serintode.c          # main(), subcommand dispatch, getopt_long
io.c / io.h          # series parser, polynomial-coeff pretty-printer (callback API), xmalloc/xcalloc
solver.c / solver.h  # kernelMP wrapper, spurious-zero check, select_best_linear, select_best_nonlin
modes_linear.c       # build_matrix_linear, print_monomial_linear, mode_descriptor mode_linear
modes_nonlin.c       # build_matrix_nonlin, print_monomial_nonlin, combs(), mode_descriptor mode_nonlin
modes_mahler.c       # (later) build_matrix_mahler, print_monomial_mahler, mode_descriptor mode_mahler
Makefile
```

`makelookup` becomes a subcommand inside `serintode.c`, sharing
`combs()` from `modes_nonlin.h`.

### Pretty-printer: callback-driven

The polynomial-coefficient block (the `(...)` with sign / zero /
leading-coeff / `x^k` logic, byte-identical across modes today) lives
once in `io.c`. Per-mode monomial output goes through a callback:

```c
typedef void (*monomial_printer)(FILE *eqs, FILE *out, long term_index,
                                 const void *mode_ctx);

void print_ode(FILE *eqs, FILE *out,
               const mpz_t *N, long numterms, long max_poly_order, long stride,
               monomial_printer print_monomial, const void *mode_ctx);
```

`numterms` is `ODE_ORDER+1` for linear/Mahler and `numterms` for nonlin.
`stride` is `nulldim`. `mode_ctx` carries mode-specific data (ODE_ORDER
for linear, `orderexp` for nonlin, `k` for Mahler).

Existing dual stdout-and-file output is preserved by passing two
`FILE*`s.

### Linear stays its own builder

Linear uses the falling-factorial trick directly on the input series
(no derivative table, no convolution). Folding it into nonlin's
builder would force allocating `D[][]` and `S[][]` and running a
singleton convolution — slower, more memory, less direct. Mahler is
shape-isomorphic to linear (one factor per term), not nonlin.

### Two named selectors

`solver.c` exposes `select_best_linear` (used by linear and Mahler) and
`select_best_nonlin` (different tiebreak chain, different acceptance
gates). The kernel call and spurious-zero check are shared.

The linear/Mahler selector minimizes `(MAX_FOUND_POLY_ORDER, termsused)`.
The nonlin selector minimizes `(finalorders, MAX_FOUND_DEPTH,
MAX_FOUND_POLY_ORDER, termsused)` and treats the
`mintermsused < finalorders+1` rejection as conditional on
`MIN_MAX_FOUND_DEPTH < 2`.

### Tunables: getopt_long now, scoped to per-run knobs

Per-run knobs become flags: `--checks`, `--min-order`, `--max-coeffs`,
`--min-depth` / `--max-depth` (nonlin only), `--lookup-dir` (nonlin
only), `--k-min` / `--k-max` (mahler only). All have defaults so a
bare `serintode <mode> <input>` works.

`MAX_LINE_LENGTH` stays a `#define` in `io.h` — it's a buffer-size
guard, not a tuning parameter.

### Output filenames: unified `<input>_<mode>_<N>-checks.txt`

`<input>_linear_<N>-checks.txt`, `<input>_nonlin_<N>-checks.txt`,
`<input>_mahler_<N>-checks.txt`. Replaces the asymmetric `_solution_`
/ `_nonlinsol_` / `_mahlersol_`. Format string lives once in `io.c`;
the suffix comes from each mode's `mode_descriptor`.

This is a one-time breaking change for any scripts globbing the old
output names.

### OOM cleanup: `xmalloc` / `xcalloc`, drop ladders

This is a one-shot computation that exits when done; the OS reclaims
everything on `exit()`. The existing nested cleanup ladders on malloc
failure (~200 LOC across the three programs) aren't doing useful work
and are subtly bug-prone (stale `fclose(foutsum)` stubs in
`serintode_iml.c`).

```c
// io.h
void *xmalloc(size_t n);            // exits on failure with strerror message
void *xcalloc(size_t n, size_t s);
```

Success-path frees stay (asymmetric, but honest).

### Backwards compat: break clean

Old binary names (`serintode_iml.o`, `serintode_iml_nonlin.o`,
`serintode_iml_nonlin_lookup.o`, `makelookup`) go away. No symlinks, no
shell stubs. README and CLAUDE.md document the new shape. Git history
preserves the old programs if anyone needs to bisect.

### Build: minimal Makefile

```make
CFLAGS  ?= -Wall -O2
IML_DIR ?= /home/ralf/math/iml
CPPFLAGS = -I$(IML_DIR)/include
LDFLAGS  = -L$(IML_DIR)/lib64
LDLIBS   = -liml -lcblas -lgmp -lm

OBJS = serintode.o io.o solver.o modes_linear.o modes_nonlin.o
# modes_mahler.o appended when that mode lands

serintode: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test: serintode
	@./serintode linear tests/central_binomials.txt --checks=6 > /dev/null
	@diff -u tests/central_binomials.txt_linear_6-checks.txt \
	         tests/expected/central_binomials_linear.txt && echo "PASS: linear"
	@./serintode nonlin tests/A037227.txt --checks=0 > /dev/null
	@diff -u tests/A037227.txt_nonlin_0-checks.txt \
	         tests/expected/A037227_nonlin.txt && echo "PASS: nonlin"

clean:
	rm -f $(OBJS) serintode

.PHONY: clean test
```

### Test harness: regression on linear and nonlin

`tests/central_binomials.txt` already has a checked-in expected output.
For nonlin, generate a baseline by running current
`serintode_iml_nonlin.o` on `tests/A037227.txt`, capture, commit. The
expected outputs live under `tests/expected/` to avoid colliding with
fresh-run output filenames.

## Migration plan (6 stages)

Each stage is one commit, leaves the tree green, and runs `make test`
before moving on.

### Stage 1: Baselines + Makefile (no code change)

- Run current `./serintode_iml.o tests/central_binomials.txt`, capture
  output, commit as `tests/expected/central_binomials_linear.txt`.
- Run current `./serintode_iml_nonlin.o tests/A037227.txt`, capture
  output, commit as `tests/expected/A037227_nonlin.txt`.
- Add Makefile that builds the *current* four programs and provides a
  `make test` target running both regressions against the old
  binaries. CLAUDE.md compile block becomes "run `make`."
- This commit catches any pre-existing flakiness before refactoring
  starts.

### Stage 2: Extract `io.c` and `solver.c`

- `io.c` / `io.h`: `xmalloc`, `xcalloc`, the series parser
  (`read_series`), the polynomial-coefficient pretty-printer
  (`print_ode` with the callback API).
- `solver.c` / `solver.h`: `kernelMP` wrapper, spurious-zero check,
  `select_best_linear`, `select_best_nonlin`.
- `serintode_iml.c` switches to the shared helpers; its OOM ladders
  drop; its monomial printing becomes a `print_monomial_linear`
  callback inline in the same file.
- The two nonlin programs are still untouched (Stage 4 cleans them up).
- `make test` passes — linear regression goes through the new shared
  code path, nonlin still uses the old binary.

### Stage 3: `modes_linear.c` + dispatcher

- Add `serintode.c` with `main()` and subcommand dispatch.
- Move the linear program's per-mode bits into `modes_linear.c`
  (matrix builder, `print_monomial_linear`, `mode_descriptor
  mode_linear`).
- `serintode linear <input>` works. Old `serintode_iml.o` Makefile
  target stays for now so the nonlin tests still have an old binary to
  compare against.
- `make test` runs `serintode linear` for the linear regression and
  the old binary for nonlin.
- Output filename is now `<input>_linear_<N>-checks.txt`. Update the
  expected-output diff target accordingly.

### Stage 4: `modes_nonlin.c` (merge nonlin + nonlin_lookup) + makelookup subcommand

- Move nonlin into `modes_nonlin.c`, collapsing the two existing
  programs. `--lookup-dir=PATH` flag toggles between in-process
  `combs()` and reading `<PATH>/o<n>d<p>.txt`. `combs()` lives in
  `modes_nonlin.c`, declared in `modes_nonlin.h`.
- Drop OOM ladders, switch to `xmalloc`/`xcalloc`.
- Wire `serintode nonlin` and `serintode makelookup` into the
  dispatcher (the latter `#include "modes_nonlin.h"` for `combs()`).
- `make test` now exercises both regressions against the new binary.

### Stage 5: getopt_long + cleanup

- Add `getopt_long` parsing in `serintode.c` for `--checks`,
  `--min-order`, `--max-coeffs`, `--min-depth`, `--max-depth`,
  `--lookup-dir`, `--k-min`, `--k-max`. Defaults match the most
  commonly-used existing values (linear: `NUM_CHECKS=6`, nonlin:
  `NUM_CHECKS=0`).
- `make test` invocations updated to pass explicit flags.
- Delete `serintode_iml.c`, `serintode_iml_nonlin.c`,
  `serintode_iml_nonlin_lookup.c`, `makelookup.c`. Remove their
  Makefile targets.
- Update `README.txt` and `CLAUDE.md` to describe the new shape only.
- Tick the Maintainability box in `TODO.md`.

### Stage 6 (separate PR): `modes_mahler.c`

- Add `modes_mahler.c` with `build_matrix_mahler`,
  `print_monomial_mahler`, `mode_descriptor mode_mahler`. Reuse
  `select_best_linear` from `solver.c`.
- Wire `serintode mahler` into the dispatcher.
- Add Mahler regression(s) to `make test` — Thue-Morse is the obvious
  candidate.
- Append `modes_mahler.o` to the Makefile's `OBJS`.

## Out of scope of this refactor

- `serintode_flint.c` (FLINT-out-of-scope per project memory).
- `old/` (archival).
- The other open `TODO.md` items (FLINT annihilation check, fscanf
  return-value check, MAX_LINE_LENGTH continuation handling, mkdir
  permission tightening) — orthogonal to the refactor and easier to
  audit one-at-a-time after the codebase has been collapsed.
