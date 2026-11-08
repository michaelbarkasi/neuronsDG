## R CMD check results

0 errors | 1 warning | 2 notes

* This is a new release (first submission to CRAN).

* The warning is a local toolchain artifact, not related to this package's
  code:

    Found the following significant warnings:
      Boolean.h:62:36: warning: unknown warning group '-Wfixed-enum-extension',
      ignored [-Wunknown-warning-option]

  This comes from Apple clang's own system header (R_ext/Boolean.h) on the
  local check machine's compiler version and does not appear on win-builder.

* One note is because `pandoc` is not installed on the local check machine
  ("Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed."); CRAN's build machines have pandoc available.

* The other note ("New submission") is expected for a first-time submission.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.2
* win-builder (R-devel, 2026-11-05 r90840 ucrt): 1 NOTE (New submission)

## Downstream dependencies

This is a new submission; there are no reverse dependencies to check.
