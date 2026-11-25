# CHO v27 research snapshot

This package publishes the wavelet-based v27 implementation used for UAB seizure research analyses. The source snapshot was checked on November 25, 2026. Those analyses explicitly set `wavelet_cycles=8` and `merge_mode='fast_tf_overlap'`, overriding the function's variable-cycle and correlation-merging defaults.

## Source and packaging

Original detector SHA-256:

```text
36a88936d140d6b09e792b4673b917462f9e0143247014e04697c2c43276f320
```

`functions/v27/CHO_v27.m` retains that detector's numerical code; its help example corrects a stale `CHO_v26` function name. Private dependencies include the original v27 aperiodic-fitting helper, wavelet v4/v6 helpers, and diagnostic line helpers. The v27 fitter is separate from the earlier v22 fitter. Helper comments are UTF-8, and temporary warning suppression restores the caller's warning state. Convolution calculations, thresholds, ACF checks, and merge rules are unchanged. Original helper attribution is retained.

`cho_v27_parameters.m` records the research preset. `example_CHO_v27.m` uses the already-public Figure 8 fixture and adds a readable two-panel visualization. UAB recordings, clinical metadata, and machine-specific paths are not distributed.

## Local validation

Fresh MATLAB R2026a Update 3 desktop sessions on macOS (Apple silicon), using the default MATLAB path plus this repository:

- Four public regression tests passed, including the legacy v22 example.
- Fixed 8-cycle v27 preset, 500 Hz, 1–150 Hz: 7 episodes on the public fixture.
- Direct v27 defaults, 500 Hz, 1–40 Hz: 5 episodes on the same fixture.
- Packaged and original local implementations returned identical complete result structures on both public configurations and one local UAB reference signal.
- The local reference had 8 episodes. Timing, frequency, and ACF fields matched its saved research result; the largest peak-value difference was approximately 2.18e-10. The saved-value comparison used a relative tolerance of 1e-9. Private input data and results remain local.
- The public result figure was opened and inspected. It displays 1–40 Hz while analyzing 1–150 Hz; all seven detections fall within the displayed range.

GitHub Actions installs the five documented toolboxes and runs the public tests on R2026a/Ubuntu. Consult the workflow result for a specific commit. Other MATLAB versions and detection performance on other recordings have not been established by this packaging check.

## Reproduce and report

Run `add_path`, then `cho_output = example_CHO_v27`. For the UAB-style configuration on another signal, pass `cho_v27_parameters()` and the actual acquisition sampling rate to `CHO_v27`.

Keep the effective `param`, actual `freqs`, acquisition sampling rate, source commit SHA, and preprocessing/channel-selection decisions. The function defaults and explicit research preset are deliberately distinct.
