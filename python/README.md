# Detect oscillatory episodes with CHO in Python

Run the detector on the single-channel hippocampal recording used in the MATLAB Figure 8 example. From the repository root, with Python 3.10 or newer:

```sh
python -m pip install "numpy>=1.22" "scipy>=1.9" "matplotlib>=3.6"
python python/cho_demo.py --out cho-demo-output
```

The script reads the bundled [raw_data signal](../data/data_for_figure8.mat) at 500 Hz. It estimates the aperiodic background across four time windows, screens elevated narrow-band envelopes for at least two cycles, checks periodicity against the original signal's autocorrelation, and merges correlated, overlapping events.

**Output.** Open `detections.png` to see the raw trace and the background-corrected time–frequency map. White boxes mark the final episodes. In `events.csv`, `start` and `stop` are zero-based inclusive sample indices (divide by 500 for seconds), `low`, `high`, and `center` are in Hz, and `peak` is the peak background-corrected log envelope.

**Scope.** This example uses the 1–40 Hz Figure 8 settings. The reference implementation is [CHO_v22.m](../functions/CHO_v22.m). The Python implementation uses SciPy's bounded linear solver for the fixed-zero-knee background fit and clips frequency neighborhoods at the edges. Event-level numerical agreement with MATLAB has not been established. See the [eLife paper](https://elifesciences.org/articles/91605) for the method.
