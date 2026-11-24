# Cyclic Homogeneous Oscillation detection (CHO)

CHO is a MATLAB method for detecting oscillatory episodes in neural recordings. It combines time-frequency analysis, aperiodic-background estimation, cycle-count criteria, and autocorrelation to identify oscillations and estimate their fundamental frequency. The [eLife publication](https://doi.org/10.7554/eLife.91605.3) describes the method and its evaluation.

The detector is implemented in [`functions/CHO_v22.m`](functions/CHO_v22.m). This repository includes a hippocampal example, sample data, and simulation scripts.

## Requirements

The hippocampal example requires MATLAB and the following toolboxes:

| Toolbox | Functions used by the example |
| --- | --- |
| Econometrics Toolbox | `autocorr` |
| Signal Processing Toolbox | `hilbert` |
| Image Processing Toolbox | `bwlabel` |
| Optimization Toolbox | `optimoptions`, `lsqcurvefit` |

The example has been tested with MATLAB R2026a Update 3 on macOS (Apple silicon), with version 26.1 of these toolboxes. Other MATLAB releases and operating systems have not been verified in this example check.

## Run the hippocampal example

Download or clone this repository, open the repository root as the current folder in MATLAB, and run:

```matlab
add_path
figure8BCDE
```

[`add_path.m`](add_path.m) adds the repository's `functions` and `data` folders to the MATLAB path. FieldTrip and EEGLAB are not required for this example. If either is already on the path, `add_path` also initializes it; scripts that use those packages require their corresponding installations.

[`figure8BCDE.m`](figure8BCDE.m) loads the included [`data/data_for_figure8.mat`](data/data_for_figure8.mat) and calls `CHO_v22` with a sampling rate of **500 Hz** and `param.frequency_axis = 1:40`. The input is a 5,000-sample, single-channel signal (10 seconds).

With plotting enabled, the detector creates a figure with four panels: the raw signal with final detections, a time-frequency map with initial candidates, a time-frequency map with final detections, and the spectrum with its aperiodic fit.

### Inspect the results

The script leaves `cho_output` in the MATLAB workspace:

```matlab
class(cho_output)
size(cho_output)
numel(cho_output.bounding_boxes)
cho_output.bounding_boxes
```

For the included example, the local check described below returned a `1-by-1 struct` with six final detected episodes. Its fields were:

| Field | Example size | Description |
| --- | --- | --- |
| `bounding_boxes` | `1-by-6 struct` | Detected episodes, with `center_fp`, `center_tp`, `peak_val`, `minF`, `maxF`, `start`, and `stop` |
| `initial_ap_fit_coef` | `1-by-3 double` | Aperiodic-model fitting coefficients |
| `init_ap_fit` | `40-by-1 double` | Aperiodic fit evaluated at the analyzed frequencies |
| `signal_std` | `1-by-1 double` | Standard deviation of the input signal |

For this example's 1 Hz frequency grid, `center_fp` gives the center frequency in Hz. `center_tp`, `start`, and `stop` are sample indices. The bounding boxes are not necessarily sorted by onset.

**Example verification (November 24, 2026).** The two commands above were run in a fresh MATLAB desktop session using code at [`532193a`](https://github.com/neurotechcenter/CHO/commit/532193a88d5203f9ce23d3508f093839978b6759). A separate session started with `restoredefaultpath` produced identical `cho_output`. The check confirmed that the detector received the included signal, 500 Hz sampling rate, and frequencies `1:40`; it returned six episodes with center frequencies between 7 and 10 Hz and created one four-panel figure. The exported figure was visually inspected. This check covers the included example, rather than all simulation scripts or a new assessment of detection accuracy.

## Simulation scripts

The simulation examples are [`figure4BCD.m`](figure4BCD.m) and [`figure4EFG.m`](figure4EFG.m).

## Cite the method

Please cite the final published article when using CHO:

Cho H, Adamek M, Willie JT, Brunner P. (2024). **Novel cyclic homogeneous oscillation detection method for high accuracy and specific characterization of neural dynamics.** *eLife*, **12**, RP91605. [https://doi.org/10.7554/eLife.91605.3](https://doi.org/10.7554/eLife.91605.3).

The Version of Record was published on September 6, 2024. The [publisher's article page](https://elifesciences.org/articles/91605) provides the version history and peer reviews; the [full text is also available in PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11379461/). The earlier [bioRxiv preprint](https://doi.org/10.1101/2023.10.04.560843) is retained as a historical reference.

```bibtex
@article{Cho2024CHO,
  author  = {Cho, Hohyun and Adamek, Markus and Willie, Jon T and Brunner, Peter},
  title   = {Novel cyclic homogeneous oscillation detection method for high accuracy and specific characterization of neural dynamics},
  journal = {eLife},
  year    = {2024},
  volume  = {12},
  pages   = {RP91605},
  doi     = {10.7554/eLife.91605.3},
  url     = {https://doi.org/10.7554/eLife.91605.3}
}
```

## Related patent

Brunner P, Cho H. **Systems and methods for detection of neurophysiological signal oscillations.** U.S. Patent **12,678,081 B2**, issued **July 14, 2026**. Assignee: **Washington University**, St. Louis, Missouri, USA. See the [USPTO Official Gazette record](https://patentsgazette.uspto.gov/week28/OG/html/1548-2/US12678081-20260714.html).

The patent was filed on March 31, 2023 as U.S. Application **18/193,930**, claiming priority to provisional application **63/326,257** filed March 31, 2022. The earlier application publication is [US 2024/0074687 A1](https://patents.google.com/patent/US20240074687A1/en), published March 7, 2024.

## License

The source code is distributed under the [MIT License](LICENSE).
