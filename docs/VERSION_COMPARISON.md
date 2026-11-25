# CHO v22–v27 comparison and technical specification

This document describes the implementations actually included in this repository. It separates **v22**, **v27 function defaults**, and the **explicit v27 UAB research preset**. The preset changes parameters; it is not another detector version. The [getting-started guide](../README.md) runs the newer research configuration on public data.

`v22` and `v27` identify MATLAB detector implementations. They are separate from repository release-tag numbering. Use the source commit SHA as well as the function name to identify an exact analysis implementation.

## Which configuration should I choose?

| Goal | Entry point | Configuration |
| --- | --- | --- |
| Reproduce the earlier Figure 8 example | `add_path; figure8BCDE` | v22, 500 Hz input, 1–40 Hz, minimum 2 cycles |
| Start with v27's built-in settings | `CHO_v27(signal, fs)` | 1–40 Hz, variable 3–10 wavelet cycles, minimum 2 cycles, correlation merging |
| Reproduce the UAB-style research settings | `CHO_v27(signal, fs, cho_v27_parameters())` | 1–150 Hz, fixed 8 wavelet cycles, minimum 3 cycles, fast merging |
| Run a complete public v27 example | `add_path; result = example_CHO_v27` | UAB-style preset, public 500 Hz fixture, result figure |

Some earlier UAB review analyses used `frequency_range=1:100`; the published preset uses `1:150`, matching the later TF/fundamental-frequency analysis. Record which range you select.

## Algorithm differences

| Component | v22 | v27 |
| --- | --- | --- |
| Spectral front end | Calls `fft_bandpass_filtering` at each requested frequency, with a band from center minus the frequency step to center plus the step; takes a Hilbert envelope | Calls the included `wavelet_single_freq_filtering_v6`; takes the magnitude of its complex Morlet convolution |
| Envelope transform | `log10(abs(hilbert(filtered_signal)))` | `log10(abs(complex_wavelet_signal))` |
| Wavelet width | Not applicable | Configurable cycles; defaults to a linear progression from 3 to 10 across the analyzed grid. The research preset explicitly selects 8 at every frequency |
| Aperiodic fit | Fits four time chunks; selects the fit with the lowest fitted value at the first analyzed frequency | Fits four chunks by default (`n_ap_fit_chunks`); selects the fit with the lowest mean fitted value across analyzed frequencies |
| Candidate threshold | Fixed `2 * std(residual_map(:))` | `osc_threshold_std * std(residual_map(:),'omitnan')`; default and preset multiplier are 2 |
| Minimum episode length | Caller supplies `minimum_cycles`; Figure 8 uses 2 | Defaults to 2; the research preset uses 3. Both implementations calculate cycles from `(stop-start)/fs` |
| ACF landmark location | Midpoints of positive and negative suprathreshold ACF regions, rounded to sample lags | Positive ACF-region maxima, with optional parabolic interpolation around each maximum; interpolation is on by default |
| ACF acceptance | Requires sufficient positive and negative landmarks; ACF frequency must lie in the candidate frequency bounds, and a positive or negative interval-consistency test must differ from 100% by less than 30 | Uses positive landmarks; defaults to frequency-bin acceptance and an interval-consistency difference from 100% of less than 10 (`conceding_threshold_for_peaks`) |
| Reported center frequency | Replaced with rounded integer-Hz positive-ACF frequency | By default replaced with the floating-point positive-ACF frequency; optional `acf_preserve_wavelet_center=1` retains the wavelet-bin center |
| Additional ACF mode | No configurable period-tolerance mode | `acf_validation_mode='period_tolerance'` is available; default tolerance is 25%. **This mode is not enabled by the function defaults or the UAB preset** |
| Merging | Requires temporal overlap and center-frequency separation of at most 2 Hz, then positive filtered-signal correlation with `p < 0.05` | Selectable correlation/overlap, fast time-frequency overlap, or no merging; see rules below |
| Input handling | No dedicated empty/flat/nonfinite-input status handling | Converts a supplied single-channel vector to a double column, interpolates NaN gaps where possible, and returns status values for empty, flat/invalid, or unavailable-frequency input |
| Frequency handling | Uses `param.frequency_axis` directly; no explicit Nyquist clipping | Uses `frequency_range`/`frequency_step`, normalizes the requested grid and clips it to a one-step-below-Nyquist guard |
| Output metadata | Final boxes, signal standard deviation, and legacy fitting outputs | Adds effective parameters/grid, status, stage counts, ACF diagnostics, selected fit, and optional maps |

Both versions use an aperiodic model of the form `b1 - log10(b2 + f^b3)` and connected components followed by cycle and autocorrelation checks. Changing the front end, fit selection, ACF rules, or merging can change detections even when the nominal frequency range and threshold multiplier match.

In the included fitters, the lower and upper bounds both fix `b2` to zero; this is not an unconstrained knee fit. The retained v27 fitter also calls `ttest2`, but its returned test statistics do not participate in the active fit/acceptance condition.

The v27 header describes support for high-frequency ACF improvements. The active defaults still use `acf_validation_mode='frequency_bin'` and `acf_preserve_wavelet_center=0`. Selecting v27 does not automatically select period-tolerance acceptance or preserve the wavelet center. Subsample ACF interpolation changes frequency quantization; it does not establish a particular frequency-estimation accuracy.

Source: [v22 front end and background selection](../functions/CHO_v22.m#L58-L103), [v22 ACF validation](../functions/CHO_v22.m#L219-L289), [v27 defaults](../functions/v27/CHO_v27.m#L67-L94), [v27 front end/background fit](../functions/v27/CHO_v27.m#L114-L164), and [v27 ACF validation](../functions/v27/CHO_v27.m#L289-L393).

### Exact merging rules

**v22:** For overlapping time intervals with center-frequency separation no greater than 2 Hz, the merge helper filters at each box's frequency bounds, applies a Tukey window, and merges when Pearson correlation is positive and its `p` value is below 0.05. Although the example supplies `ovlp_threshold=0.5`, that variable is not used by the active v22 merge-acceptance condition.

**v27 correlation mode:** The helper uses wavelet-v4 real filtered signals and a Tukey window. After requiring temporal overlap, it merges if either positive correlation has `p < 0.05`, or the time-frequency intersection exceeds `ovlp_threshold` times either box's volume. The v22 2-Hz center-separation rejection is commented out in this v27 helper. If more than `max_boxes_for_correlation_merge` candidates reach merging (default 250), v27 takes the fast path even when correlation mode was requested.

**v27 fast mode:** Two boxes must overlap in time. Their inclusive time overlap, divided by the shorter box's inclusive duration, must be at least `fast_merge_time_overlap_frac`. They must also satisfy the implemented frequency-gap or center-proximity test. Defaults/preset use a time fraction of 0.05 and frequency-gap setting of 1 Hz. The frequency test uses `max(minF) - min(maxF)` or center separation no greater than `max(fast_merge_freq_gap_hz, frequency_step)`. Connected merge groups are combined. `ovlp_threshold` does not control this fast path.

Source: [v22 merge helper](../functions/CHO_v22.m#L500-L550), [v27 mode selection](../functions/v27/CHO_v27.m#L441-L490), [fast predicate](../functions/v27/CHO_v27.m#L686-L700), and [v27 correlation helper](../functions/v27/CHO_v27.m#L768-L817).

## Input and numerical specification

| Item | Specification |
| --- | --- |
| Signal | A real, single-channel time series; use a time-by-one numeric vector. Select a channel before calling the detector. v27's `signal(:)` conversion will flatten a matrix rather than perform a multichannel analysis |
| Sampling rate | Caller-supplied acquisition sampling rate in Hz. The algorithm is not restricted to 500 Hz; the public example uses 500 Hz. Supply a positive, finite value |
| Frequency grid | v22 example: `frequency_axis=1:40`. v27 defaults: `frequency_range=1:40`, `frequency_step=1`. Research preset: `1:150`, step 1 Hz |
| v27 frequency guard | Frequencies are limited to `floor(fs/2)-frequency_step`; for the 500 Hz, 1-Hz-step example the upper guard is 249 Hz. A requested 1–150 Hz range is clipped when the sampling rate is too low |
| Arbitrary v22 axes | Legacy code uses frequency-array indices as Hz in parts of the detector. The included consecutive 1-Hz grid is tested; arbitrary non-unit/nonconsecutive grids are not established by that test |
| Wavelet support | Width scales with cycles divided by frequency. Low-frequency wavelets therefore span longer time intervals; the included helpers use FFT convolution. Edge effects and window length matter when interpreting detections |
| Frequency precision | Grid spacing is a configuration value, not a measured accuracy guarantee. Report the grid, sampling rate, wavelet cycles and ACF settings together |
| Time convention | Box `start`, `stop`, and `center_tp` are one-based sample indices. Onset in seconds is `(start-1)/fs`; endpoint separation is `(stop-start)/fs`, while inclusive sample count is `stop-start+1` |
| Magnitude convention | Candidate and `peak_val` values use residual log10 envelope/wavelet magnitude after the fitted aperiodic background is subtracted. These are not raw signal amplitudes or dB power values |
| Nonfinite v27 input | NaN gaps are linearly filled with nearest endpoint values where possible. If nonfinite values remain, or the signal is flat, the detector returns an invalid/flat-input status. Review input quality before interpretation |
| Record length | No universal minimum duration is specified or validated here. Cycle requirements and wavelet support depend on frequency. These tests do not establish behavior for every short recording or parameter combination |

The bundled example receives the public fixture directly; it does not add resampling, rereferencing, or artifact rejection. Apply and report any preprocessing required by your own research pipeline.

## Parameter reference

These are effective settings from the source, not suggested replacements for the research configuration.

| Parameter | v22 Figure 8 example / implementation | v27 default | v27 research preset |
| --- | --- | --- | --- |
| `plot` | Example sets 1; caller supplies it | 0 | 0; public wrapper makes its own figure |
| Frequency argument | `frequency_axis=1:40` | `frequency_range=1:40` | `frequency_range=1:150` |
| `frequency_step` | Inferred from first two axis entries | 1 Hz | 1 Hz |
| `wavelet_cycles` | Not used | `'linear_3_10'` | 8 |
| `minimum_cycles` | Example sets 2 | 2 | 3 |
| `osc_threshold_std` | Not configurable; multiplier 2 | 2 | 2 |
| `n_ap_fit_chunks` | Fixed 4 | 4 | 4, inherited |
| `AC_PK_Threshold` | Forced to 1 in the function (`autocorr` NumSTD) | 1 | 1, inherited |
| `conceding_threshold_for_peaks` | Fixed interval check of less than 30 percentage points | 10 percentage points | 10, inherited |
| `acf_validation_mode` | Frequency-bin test | `'frequency_bin'` | `'frequency_bin'`, inherited |
| `acf_period_tolerance_pct` | Not used | 25; active only in period-tolerance mode | Inherited but inactive in the selected mode |
| `acf_lag_interpolation` | Not configurable; rounded region midpoint | 1 | 1, inherited |
| `acf_preserve_wavelet_center` | Not used | 0 | 0, inherited |
| `merge_mode` | Fixed correlation helper | `'correlation'` | `'fast_tf_overlap'` |
| `max_boxes_for_correlation_merge` | Not used | 250 | Inherited; explicit fast mode is already selected |
| `ovlp_threshold` | Example passes 0.5; inactive in active v22 acceptance condition | 0.5; used in correlation/volume-overlap mode | 0.5, retained but inactive in fast mode |
| `fast_merge_time_overlap_frac` | Not used | 0.05 | 0.05 |
| `fast_merge_freq_gap_hz` | Not used | `max(1, frequency_step)` | 1 Hz |
| `return_maps` | No equivalent output switch | 0 | Inherited 0; public wrapper sets 1 for visualization |
| `birdsong_tfmap`, `plot_birdsong` | Not used | 0, 0 | 0, 0 |
| `three_d_plot`, `acf_plots` | Forced to 0 in the function | 0, 0 | 0, 0, inherited |

Optional birdsong-style display defaults are `birdsong_freq_range=[]` (use analyzed endpoints), `birdsong_n_freq=80`, `birdsong_baseline_idx=[]` (use the full window), `birdsong_clim=[-6 12]`, and `birdsong_colormap='turbo'`. They control that optional display map; the research preset leaves this path disabled.

In v27, missing or empty parameter fields take their defaults. The research preset explicitly sets only selected fields; other values are inherited from the detector. Save the effective `result.param`, including `wavelet_cycles_by_freq`, rather than only the smaller input preset.

## Output specification

Both versions return a scalar result structure for the supplied channel/window.

| Output | v22 | v27 |
| --- | --- | --- |
| `bounding_boxes` | Final episode array | Final episode array; empty input/status cases can contain an empty numeric array, so check `isempty` before indexing |
| Common box fields | `start`, `stop`, `center_tp`, `center_fp`, `minF`, `maxF`, `peak_val` | Same fields |
| ACF box diagnostics | Not included | `acf_fp`, `acf_period_error_pct`, `acf_interval_similarity_pct` |
| `signal_std` | Input signal standard deviation | Standard deviation after the detector's input conversion/gap handling; may be NaN for early status returns |
| Legacy fitting fields | `initial_ap_fit_coef`, `init_ap_fit` | Retained for compatibility |
| Selected background | Not returned in a dedicated `ap_fit` field | `ap_fit` is the selected background actually subtracted |
| Status and settings | No equivalent status/grid/parameter fields | `status`, `freqs`, `param` |
| Stage counts | Not returned as named counts | `n_initial_boxes`, `n_after_cycle_check_boxes`, `n_before_autocorr_check_boxes`, `n_after_autocorr_check_boxes`, `n_final_boxes` |
| Intermediate boxes | Not returned | `initial_bounding_boxes`, `bounding_boxes_second_step`, `before_auto_corr_check` |
| Optional maps | Not returned | `filtered_data`, `env_data`, `denoised_env_data`, `masked_env_data`, `oevent_map` when `return_maps=1`; arrays are time-by-frequency |
| Optional display TF map | Not returned | `birdsong_tfmap` and optional figure handle when requested; this display path uses its own 3–10-cycle map and does not replace the detector's configured wavelet cycles |

The legacy fitting fields are assigned from the last chunk's fitter output. They should not be assumed to describe the selected background chunk; in v27 use `ap_fit` for the selected curve. Successful normal v27 results do not include a guaranteed `srate` field; save the acquisition sampling rate separately. Status values include `"ok"`, `"empty_signal"`, `"invalid_or_flat_signal"`, and `"empty_frequency_range_after_nyquist_clip"`.

## MATLAB and dependency specification

Local verification used MATLAB R2026a Update 3 on Apple silicon; GitHub Actions uses R2026a on Ubuntu. Other releases have not been established by these checks.

| Toolbox | v22 use | v27 use |
| --- | --- | --- |
| Econometrics | `autocorr` | `autocorr` |
| Signal Processing | `hilbert`, `tukeywin` | `tukeywin` in correlation merging |
| Image Processing | `bwlabel` | `bwlabel` |
| Optimization | `optimoptions`, `lsqcurvefit` | `optimoptions`, `lsqcurvefit` |
| Statistics and Machine Learning | `corr` in merging | `corr` in correlation merging and `ttest2` in the retained v27 fitter |

The five-toolbox package is the documented/tested installation. Wavelet Toolbox is not required by these custom FFT-based Morlet helpers. Public examples run without FieldTrip, EEGLAB, or Parallel Computing Toolbox. The detector call itself is single-channel; batching and any parallel workers belong to the surrounding analysis script. Hardware limits, throughput, GPU support, and maximum record size have not been benchmarked for this package.

## Migration examples

Reproduce v22 with its own argument name:

```matlab
add_path
p22 = struct('plot', 0, 'minimum_cycles', 2, ...
    'ovlp_threshold', 0.5, 'frequency_axis', 1:40);
r22 = CHO_v22(signal, fs, p22);
```

Use the explicit v27 research preset:

```matlab
p27 = cho_v27_parameters();
p27.frequency_range = 1:150;
r27 = CHO_v27(signal, fs, p27);
save('cho_result.mat', 'r27', 'fs', 'p27')
```

For a controlled frequency-range comparison, set `p27.frequency_range=1:40`. That aligns the requested range, while the front end, ACF rules, merge mode, and minimum-cycle setting still differ. Do not obtain v27 settings merely by passing a v22 `param` structure: `frequency_axis` does not configure v27's `frequency_range`.

## What has been verified?

| Public 10-second hippocampal fixture, 500 Hz | Analysis settings | Final episode count |
| --- | --- | --- |
| v22 example | 1–40 Hz; minimum 2 cycles | 6 |
| v27 direct defaults | 1–40 Hz; variable 3–10 wavelet cycles; minimum 2 cycles | 5 |
| v27 research-preset example | 1–150 Hz; fixed 8 wavelet cycles; minimum 3 cycles; fast merging | 7 |

Four regression tests cover the legacy example, v27 preset and defaults, figure/output contracts, warning-state preservation, and flat-input/Nyquist behavior. The detector/fixture revision was verified at [commit 95da794](https://github.com/neurotechcenter/CHO/commit/95da7941c038ce962722bc8681bf9b0ece90e863), including [Linux CI](https://github.com/neurotechcenter/CHO/actions/runs/45651296349). See [source provenance](CHO_V27.md) for the local research-reference comparison.

These counts come from different configurations on a single fixture. They are not a matched sensitivity/specificity comparison, a speed benchmark, or evidence that one version performs better on every recording. Share the version, source commit, effective settings, sampling rate, and preprocessing with any reported result.
