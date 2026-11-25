# Contributing to CHO

Use a separate branch and submit a pull request against `main`. Keep each change focused so reviewers can assess the implementation and its effect on research results.

## Set up and test

Install MATLAB with Econometrics Toolbox, Signal Processing Toolbox, Image Processing Toolbox, Optimization Toolbox, and Statistics and Machine Learning Toolbox. From the repository root in a fresh MATLAB session, run:

```matlab
results = runtests('tests');
assertSuccess(results);
```

Tests run the v27 research preset and function defaults, check edge cases and warning-state preservation, and retain the v22 `figure8BCDE` regression. They check outputs and figures and write PNG/MAT artifacts to `artifacts/`. Cleanup restores session settings and closes only figures created by the tests. The v22 script's `clear` command is isolated from the caller's workspace.

GitHub Actions runs these tests on Ubuntu with MATLAB R2026a and the required toolboxes. These checks cover examples and output contracts; they do not reproduce the full paper's simulation suite or establish detection performance. Keep private recordings and identifiers out of public fixtures and logs.

## Report a problem

Include the commit SHA, operating system, MATLAB release, relevant toolbox versions, exact commands, complete error text, and expected versus actual behavior. For plotting problems, include a figure or screenshot when possible.

## Propose a change

Explain the problem and resulting behavior in the pull request. State exactly what you tested and distinguish a local MATLAB run from a completed CI run.

For detector, threshold, or data changes, explain the scientific rationale and provide before/after results. Update the regression expectations only when the changed behavior is intentional and supported by evidence. A dependency or documentation fix should not require changing the detector to make its test pass.

Preserve existing data files and historical analysis scripts unless their modification is part of the proposed change. Do not commit generated simulation MAT files, personal MATLAB paths, or credentials.

## Citation and licensing

Keep the preferred article citation in `CITATION.cff` consistent with the README. Source-code contributions use the repository's existing MIT License. Do not add a version number or tag until a release is intentionally prepared.
