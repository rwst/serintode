# TODO

Findings from an audit of the C sources. None are exploitable through the
data file — they are robustness, correctness, and quality items.

## Correctness

- [ ] `serintode_flint.c` does not verify that the returned null-space vectors
      actually annihilate the input matrix. README and the file's own header
      comment note that this is why the FLINT version "can give spurious
      results unless enough checks are given." Port IML's annihilation check
      (`kernelMP` / `nullspaceMP` already does it internally) into the FLINT
      path before treating a null vector as a candidate ODE.

## Robustness

- [x] Double-free / use-after-close in the "couldn't open output equations
      file" error path: removed the redundant `fclose(fin)` / `free(M)` /
      `mpz_clear` cleanup from the error branch (the `exit(EXIT_FAILURE)`
      reclaims everything anyway). `-Wuse-after-free` warnings now silent in
      all three IML programs.

- [x] Check return value of `mpz_init_set_str` / `mpz_set_str` when parsing
      each line. Errors out with the offending line number; previously a
      malformed line silently became a zero coefficient.

- [x] Check `fscanf` return value when reading lookup tables; a truncated file
      previously left `orderexp[i][j]` undefined. Now errors out at the
      offending `(i, j)` with the lookup file path. See `modes_nonlin.c:84`
      (in `load_orderexp_from_lookup`).

- [ ] Detect lines longer than `MAX_LINE_LENGTH-1` (no newline read by
      `fgets`) and either error out or join the continuation. Currently such a
      line is silently split into two coefficients. Only matters for
      coefficients with ≥99999 digits.

- [x] Moved `input_string` from a `MAX_LINE_LENGTH+1` stack array to a
      heap buffer (free'd right after `fclose(fin)`). Now lives once in
      `read_series` (`io.c`). FLINT version out of scope.

## Latent buffer overflows (not reachable from input data)

- [x] Replace `sprintf` into the 64-byte `fouteqsname` / `fname` stack buffers
      with `snprintf` in the in-scope IML files and `makelookup.c`. (The
      `serintode_flint.c` site is out of scope — see project memory.)

## `makelookup.c`

- [x] Validate `atol(argv[1])` / `atol(argv[2])` — rejects ≤0 with a clear
      error.
- [x] `(dcheck==-1)&(errno!=EEXIST)` → `&&` at `makelookup.c:72,81`.
- [ ] `mkdir(dirname, S_IRWXU|S_IRWXG|S_IRWXO)` creates the lookup dir as
      0777. Tighten to 0755 (or 0700) unless world-writable is intentional.
  - `modes_makelookup.c:33`

## Ergonomics

- [x] Take `finname` from positional `argv[1]` and per-run knobs
      (`--checks`, `--min-order`, `--max-coeffs`, `--min-depth`,
      `--max-depth`, `--lookup-dir`, `--k-min`, `--k-max`, `--k`,
      `--max-depth`) from `getopt_long`. Bare invocation prints a
      mode-specific usage line. `MAX_LINE_LENGTH` is the only structural
      limit that still needs a recompile (`#define` in `io.h`).

## Maintainability

- [x] The four `serintode_*` programs each carry their own copy of the input
      reader, the output writer, and the main null-space loop. Refactored
      into a single dispatcher binary `serintode` with subcommands `linear`,
      `nonlin`, `makelookup`, sharing `io.{c,h}` (parser, pretty-printer)
      and `solver.{c,h}` (kernel wrapper, null-vector selectors). Per-run
      knobs are `getopt_long` flags. See `plan.md` for the design notes.

## Mahler / k-regular search — implemented, follow-ups remain

Implemented as `modes_mahler.c` (the `serintode mahler` subcommand) and
`modes_kkernel.c` (the `serintode kkernel` diagnostic). Captures
k-regular sequences (Allouche-Shallit) and, with the extra finite-image
condition reported by the `*Image size: m of N terms*` banner line,
k-automatic sequences. Out of reach of the linear / nonlin programs:
automatic sequences are not D-finite over `Q` unless eventually periodic
(Adamczewski-Bell).

The pure-Python prototype `mahler_guess.py` and its caveats in
`TODO_mahler.md` were the original sketch.

### Math (background, still useful as documentation)

Given `a_0, ..., a_{N-1}`, look for `p_0, ..., p_d in Z[x]` of degree
`<= D`, not all zero, with

    sum_{i=0}^{d} p_i(x) f(x^(k^i)) = 0,    f(x) = sum a_n x^n.

Row `N` of the matrix has entry `M[N][(i,j)] = a_{(N-j)/k^i}` when
`k^i | (N-j)` and the index is in range, else 0. Bound formulas:

    MAX_POLY_ORDER(n)  = floor((NUM_COEFFS - NUM_CHECKS) / (n+1)) - 1
    MAX_MAHLER_ORDER   = floor(log_k(NUM_COEFFS - NUM_CHECKS))

The IML precision constraint
`ceil(rows/2)*(p-1)^2 + (p-1) <= 2^53 - 1` carries over unchanged.

### Follow-ups

- [ ] **Inhomogeneous Mahler equations.** The homogeneous form misses
      the natural `{0,1}`-Thue-Morse, which satisfies
      `f(x) - (1-x) f(x^2) = x/(1-x^2)`. Either augment the matrix with
      columns for a small basis of rational inhomogeneities
      (truncations of `1/(1-x)^t`, `x^a/(1-x^b)` for small a,b,t), or
      solve for the difference of two augmented-system solutions.

- [x] **k-kernel rank diagnostic.** Standalone subcommand
      `serintode kkernel <input> [--k=N] [--max-depth=N]` that prints
      the rank of the matrix with rows `[a(k^i n + j)]_n` for each
      depth, so the user can see whether the rank stabilises (k-regular)
      or keeps growing (not k-regular). Implemented in `modes_kkernel.c`
      via `compute_rank` (cols − nullity from IML's `kernelMP`).

- [ ] **Auto pre-filter the Mahler search.** Wire `kkernel` into the
      `mahler` subcommand so each candidate `k` is rejected up-front if
      its rank keeps growing within the data's resolution. Only worth
      adding if the Mahler search becomes the bottleneck.

- [ ] **`K(x)`-coefficient lift.** Clear denominators to extend the
      search to k-Mahler equations with rational-function coefficients.
      Degree blows up; only do this if the integer-coefficient sweep
      finds nothing.

- [ ] **Cross-check vs. `mahler_guess.py`.** Diff C hits against the
      Python prototype on a shared OEIS slice to catch implementation
      bugs in either. Note: TODO_mahler.md's claim that Stern's k=2
      kernel rank is 5 disagrees with both the C `kkernel` (rank 2) and
      direct hand-calculation; one of TODO_mahler / Allouche-Shallit /
      our implementation is wrong, worth resolving as part of the
      cross-check.

## Out of scope of this audit

- [ ] `old/serintode_iml.c`, `old/serintode_nonlin2.c`,
      `old/serintode_nonlin_twoterms.c` were not audited. They likely share
      the same parser-return and `sprintf`-into-fixed-buffer issues. Either
      audit them or delete the directory if it's truly archival.
