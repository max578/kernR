# Contributing to kernR

Thank you for your interest in improving kernR. Contributions of all
kinds are welcome: bug reports, feature requests, documentation, and
code.

## Reporting issues

Open an issue at <https://github.com/max578/kernR/issues>. For bugs,
please include a minimal reproducible example, the output of
`sessionInfo()`, and the kernR version.

## Pull requests

1. Fork the repository and create a feature branch from `main`
   (`feature/<short-description>` or `bugfix/<short-description>`).
2. Install development dependencies and confirm a clean check:
   ```r
   devtools::install_dev_deps()
   devtools::check()   # must pass with no errors or warnings
   ```
3. Add or update tests under `tests/testthat/` for any behavioural
   change. Stochastic tests must set a `seed`.
4. Add a `NEWS.md` bullet for any user-facing change.
5. Keep the diff focused; do not restyle unrelated code.
6. Reference the issue your PR closes in the PR body, e.g.
   `Fixes #12`.

## Code style

Code follows a hand-crafted, reviewer-readable R style: `snake_case`
names, `<-` assignment, `TRUE`/`FALSE`, two-space indentation, an
80-character soft margin, and `pkg::fun` qualification for imported
functions. The repository carries an `air.toml`; run the
[Air](https://posit-dev.github.io/air/) formatter on code you touch.
Roxygen2 (markdown) documents every exported function.

## Code of conduct

By participating you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).
