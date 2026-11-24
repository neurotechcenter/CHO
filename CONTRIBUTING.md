# Contributing to CHO

Use a separate branch and submit a pull request against `main`. Keep each change focused so reviewers can assess the implementation and its effect on research results.

## Set up and test

Install MATLAB with Econometrics Toolbox, Signal Processing Toolbox, Image Processing Toolbox, and Optimization Toolbox. From the repository root in a fresh MATLAB session, run:

```matlab
results = runtests('tests');
assertSuccess(results);
```

The regression test executes `add_path` and `figure8BCDE` in an isolated function workspace, checks the bundled signal and output, verifies the new figure, and writes a PNG and MAT result to `artifacts/`. It restores the current folder, search path, default figure visibility, and closes only figures created by the test. The public example's `clear` command does not clear the caller's workspace.

GitHub Actions runs the same test on Ubuntu with MATLAB R2026a and the required toolboxes, and retains test results and the figure as downloadable artifacts. This check covers the included hippocampal example; it does not reproduce the full paper's simulation suite.

## Report a problem

Include the commit SHA, operating system, MATLAB release, relevant toolbox versions, exact commands, complete error text, and expected versus actual behavior. For plotting problems, include a figure or screenshot when possible.

## Propose a change

Explain the problem and resulting behavior in the pull request. State exactly what you tested and distinguish a local MATLAB run from a completed CI run.

For detector, threshold, or data changes, explain the scientific rationale and provide before/after results. Update the regression expectations only when the changed behavior is intentional and supported by evidence. A dependency or documentation fix should not require changing the detector to make its test pass.

Preserve existing data files and historical analysis scripts unless their modification is part of the proposed change. Do not commit generated simulation MAT files, personal MATLAB paths, or credentials.

## Citation and licensing

Keep the preferred article citation in `CITATION.cff` consistent with the README. Source-code contributions use the repository's existing MIT License. Do not add a version number or tag until a release is intentionally prepared.
