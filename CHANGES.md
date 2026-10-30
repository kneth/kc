# Changes

## Unreleased

The 1990s sources have been modernized so kc builds cleanly with
current C compilers, and an automated test suite has been added.

### Build and tooling

* CMake replaces `kc-inst` and the per-platform Makefiles. Only a C
  compiler, flex and bison are needed; the parser and lexer are
  generated at build time.
* The sources build as C17 (or C23) without warnings under
  `-Wall -Wextra -Wpedantic` on current GCC and clang.
* `kci` runs kc in mode 3, compiles the model and runs the simulation
  in a temporary directory, so it no longer overwrites `model.c` and
  `model.h` in the current directory.
* C99 complex numbers replace the custom `Complex` type. `complex.h` is
  no longer installed.
* Floating-point contraction is disabled, so mode 1 results no longer
  depend on the compiler or CPU.
* Legacy material (retired generators, DOS copies, platform Makefiles,
  site-specific scripts) has been removed. It remains in git history.

### Tests

* A golden-output test suite runs every example model in all five
  modes and compares all generated files, including the end state of
  mode 3 simulations. Models that kc must reject are tested too.
* GitHub Actions runs the suite on Linux (GCC, clang, C23) and macOS
  (clang, GCC 14).

### Changed results

* Mode 1 uses the intended default precision (`epsa` = 1.0e-18) when
  the model sets none. Before, the precision was 0: eight models
  reported a failure in `hqr2`, and D3 printed all-zero eigenvectors.
  The eigenvector output (mostly normalization and sign) changes for
  16 models.

### Bug fixes

* kc crashed with a segfault or "double free or corruption" on many
  models with glibc 2.43, because new expression tree nodes were not
  initialized.
* An unclosed comment made kc loop forever, and a comment ending in
  `**/` was not always recognized as closed.
* Line numbers in kc's messages were too low after multi-line
  comments.
* Names longer than 24 characters overflowed a buffer; they are now
  reported as an error.
* Mode 1: fixed out-of-bounds reads when normalizing eigenvectors and
  in the `hqr2` loop exit test (`m==1` instead of `m==l`), which made
  eigenvectors depend on heap contents.
* Mode 2: the generated `reac_F` is declared as returning `double`; it
  used to truncate its result to `int`.
* Mode 3: the species names of non-autonomous models were stored one
  slot too low. Generated simulators now compile as C17.
* Mode 5: the range check of preference values could never trigger.
* Mode 5: integral constants too large for an `int` (e.g. 4e11) were
  written as 2147483647 or -2147483648, depending on the CPU. They are
  now written in exponent notation.
* Several uninitialized variables have been fixed, including the
  result of `TreeEval` and table lookups of unknown names.
* The warning for an unknown command-line option printed the wrong
  character.
