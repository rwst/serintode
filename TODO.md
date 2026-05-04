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

- [ ] Check `fscanf` return value when reading lookup tables; a truncated file
      leaves `orderexp[i][j]` undefined.
  - `serintode_iml_nonlin_lookup.c:320`

- [ ] Detect lines longer than `MAX_LINE_LENGTH-1` (no newline read by
      `fgets`) and either error out or join the continuation. Currently such a
      line is silently split into two coefficients. Only matters for
      coefficients with ≥99999 digits.

- [x] Moved `input_string` from a `MAX_LINE_LENGTH+1` stack array to a
      `malloc`'d buffer (free'd right after `fclose(fin)`) in the three IML
      programs. (FLINT version out of scope.)

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
  - `makelookup.c:71`

## Ergonomics

- [x] Take `finname` from `argv[1]` rather than recompiling. Done in all
      three IML programs; bare invocation now prints
      `Usage: <prog> <input-file>` and exits non-zero. Other tunables
      (`NUM_CHECKS`, `MIN_ODE_ORDER`, `MAX_COEFFS`, depth bounds) still
      require recompiling — escalate to `getopt` if that becomes painful.

## Maintainability

- [x] The four `serintode_*` programs each carry their own copy of the input
      reader, the output writer, and the main null-space loop. Refactored
      into a single dispatcher binary `serintode` with subcommands `linear`,
      `nonlin`, `makelookup`, sharing `io.{c,h}` (parser, pretty-printer)
      and `solver.{c,h}` (kernel wrapper, null-vector selectors). Per-run
      knobs are `getopt_long` flags. See `plan.md` for the design notes.

## New module: `serintode_iml_mahler.c` (Mahler / k-regular search)

Port the matrix-nullspace template from `serintode_iml.c` to the Mahler
operator `M_k : f(x) -> f(x^k)`, the structural detector for k-regular
sequences (Allouche-Shallit). Captures k-automatic sequences when the
additional finite-image constraint holds. Out of reach of the existing
linear / nonlin programs: automatic sequences are not D-finite over `Q`
unless eventually periodic (Adamczewski-Bell), so `serintode_iml.c`
returns nothing on Thue-Morse, Stern, regular paperfolding, etc.

A pure-Python prototype with a working Thue-Morse demo lives at
`mahler_guess.py` in this repo; caveats for that sketch are in
`TODO_mahler.md`.

### Matrix construction

Given `a_0, ..., a_{N-1}`, look for `p_0, ..., p_d in Z[x]` of degree
`<= D`, not all zero, with

    sum_{i=0}^{d} p_i(x) f(x^(k^i)) = 0,    f(x) = sum a_n x^n.

Unknowns `c_{i,j}` for `0 <= i <= d`, `0 <= j <= D` — `(d+1)(D+1)` total.
Row `N` (coefficient of `x^N` in the LHS) gives one equation:

    sum_{i,j} c_{i,j} * [k^i | (N-j)] * a_{(N-j)/k^i} = 0.

So `M[N][(i,j)] = a_{(N-j)/k^i}` when `k^i | (N-j)` and the index is in
range, else 0. Pure integer arithmetic — feeds straight into `kernelMP`
/ `nullspaceMP`. No derivatives, no `combs()`; simpler than either
existing IML program.

### Search loop

Outer loop over base `k`. By Cobham, a sequence can be k-automatic for
at most one multiplicative class, so just take the first hit:

```c
for (long k = K_MIN; k <= K_MAX; k++)
    for (long n = MIN_MAHLER_ORDER; n <= MAX_MAHLER_ORDER; n++)
        for (long deg = 0; deg <= MAX_POLY_ORDER(n); deg++)
            // build M, nullspace, pick best null vector by
            // (max poly degree, then nnz) as in serintode_iml.c
```

Bound formulas, parallel to README's derivations for the linear case:

    MAX_POLY_ORDER(n)  = floor((NUM_COEFFS - NUM_CHECKS) / (n+1)) - 1
    MAX_MAHLER_ORDER   = floor(log_k(NUM_COEFFS - NUM_CHECKS))

The IML precision constraint
`ceil(rows/2)*(p-1)^2 + (p-1) <= 2^53 - 1` carries over unchanged — IML
doesn't care that the matrix came from a Mahler operator.

### Output

Maple syntax, joining `(<poly>)*f(x^(k^i))` with ` + `. Example for
Thue-Morse:

    (1)*f(x) + (-1+x)*f(x^2) = 0

Filename suffix `_mahlersol_<NUM_CHECKS>-checks.txt`, alongside stdout
print, matching the existing programs.

### Verification

- `nullspaceMP` already checks the null vector against the matrix —
  reuse, same as the linear program.
- Plus `NUM_CHECKS` over-determination: build with
  `NUM_COEFFS - NUM_CHECKS` rows, evaluate the recovered relation
  against the held-out tail, reject on mismatch. Same pattern as
  `serintode_iml.c`.

### CLI

    ./serintode_iml_mahler.o <input-file> [k_min] [k_max]

defaults `k_min=2, k_max=10`. `NUM_CHECKS`, `MIN_MAHLER_ORDER`,
`MAX_COEFFS` stay `const` locals at the top of `main()` until the
`getopt` escalation under Ergonomics happens.

### Follow-ups

- [ ] **Inhomogeneous Mahler equations.** The homogeneous form misses
      the natural `{0,1}`-Thue-Morse, which satisfies
      `f(x) - (1-x) f(x^2) = x/(1-x^2)`. Either augment the matrix with
      columns for a small basis of rational inhomogeneities
      (truncations of `1/(1-x)^t`, `x^a/(1-x^b)` for small a,b,t), or
      solve for the difference of two augmented-system solutions.

- [ ] **k-kernel pre-filter.** Compute the rank of the matrix with rows
      `[a(k^i n + j)]_n` for `0 <= i <= K_DEPTH`, `0 <= j < k^i`. If the
      rank keeps growing as `K_DEPTH` does, the sequence is provably not
      k-regular within the data's resolution — skip the Mahler search
      for that `k`. IML computes this rank directly. Only worth adding
      if the Mahler search becomes the bottleneck.

- [ ] **`K(x)`-coefficient lift.** Clear denominators to extend the
      search to k-Mahler equations with rational-function coefficients.
      Degree blows up; only do this if the integer-coefficient sweep
      finds nothing.

- [ ] **Cross-check vs. `mahler_guess.py`.** Once the C version works,
      diff its hits against the Python prototype on a shared OEIS slice
      to catch implementation bugs in either.

### Where to copy from

`serintode_iml.c` is the right template — same shape (linear, no
`combs()`), same null-vector selection rule, same I/O. Expect comparable
LOC. The Maintainability note above applies: this would be the fourth
program with its own copy of the file parser and output writer; weigh
against the long-standing plan to factor those out before adding a
fourth offender.

## Out of scope of this audit

- [ ] `old/serintode_iml.c`, `old/serintode_nonlin2.c`,
      `old/serintode_nonlin_twoterms.c` were not audited. They likely share
      the same parser-return and `sprintf`-into-fixed-buffer issues. Either
      audit them or delete the directory if it's truly archival.
