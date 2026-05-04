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

- [ ] Take `finname`, `NUM_CHECKS`, `MIN_ODE_ORDER`, `MAX_COEFFS` from `argv`
      / a config rather than recompiling. README already calls this out as a
      known limitation.

## Maintainability

- [ ] The four `serintode_*` programs each carry their own copy of the input
      reader, the output writer, and the main null-space loop (~500–900 LOC
      each). Every fix in this file has to be applied 3–4 times. Factor the
      shared pieces (file parsing, ODE pretty-printing, IML/FLINT wrappers)
      into a small library or a single program with a backend switch.

## Out of scope of this audit

- [ ] `old/serintode_iml.c`, `old/serintode_nonlin2.c`,
      `old/serintode_nonlin_twoterms.c` were not audited. They likely share
      the same parser-return and `sprintf`-into-fixed-buffer issues. Either
      audit them or delete the directory if it's truly archival.
