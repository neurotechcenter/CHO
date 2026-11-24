"""Run the Figure 8 single-channel CHO example without MATLAB."""
import argparse
import csv
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.io import loadmat
from scipy.optimize import lsq_linear
from scipy.stats import pearsonr


@dataclass
class Box:
    start: int  # zero-based, inclusive sample indices
    stop: int
    low: int
    high: int
    center: int
    peak: float


def runs(mask):
    edges = np.flatnonzero(np.diff(np.r_[False, mask, False]))
    return [(start, end - 1) for start, end in edges.reshape(-1, 2)]


def background(freqs, spectrum):
    # With the knee fixed to zero in oof_fitting_v2, the fit is linear in log10(f).
    design = np.column_stack([np.ones(len(freqs)), -np.log10(freqs)])
    bounds = ([spectrum.min() - 2, 0], [spectrum.max() + 2, 1000])
    fit = lsq_linear(design, spectrum, bounds=bounds).x
    residual = np.maximum(spectrum - design @ fit, 0)
    keep = residual <= np.percentile(residual, 10, method="hazen")
    if keep.sum() > 3:
        return design @ lsq_linear(design[keep], spectrum[keep], bounds=bounds).x


def confirmed_frequency(segment, fs, freq):
    centered = segment - segment.mean()
    energy = centered @ centered
    if energy == 0:
        return None
    acf = signal.correlate(centered, centered, method="fft")[len(segment)-1:] / energy
    bound = 1 / np.sqrt(len(segment))
    positive = [(a+b+1)//2 for a,b in runs(acf > bound)][1:]
    negative = [(a+b+1)//2 for a,b in runs(acf < -bound)]
    if len(positive) < 2 or len(negative) < 2:
        return None
    period = positive[0]
    hz = int(np.floor(fs / period + 0.5))
    regular = any(abs((lags[1]-lags[0])/period - 1) < 0.30 for lags in (positive, negative))
    if regular and 1 <= hz <= 40 and abs(hz - freq) <= 1:
        return hz


def detect(raw, fs=500):
    raw = np.asarray(raw, dtype=float).squeeze()
    if raw.ndim != 1 or len(raw) < fs or not np.isfinite(raw).all() or np.ptp(raw) == 0:
        raise ValueError("Expected at least one second of finite, varying single-channel samples.")
    freqs = np.arange(1, 41)
    n = len(raw)
    fft = np.fft.fft(raw)
    bins = np.arange(n//2 + 1) * fs / n

    def band(low, high):
        # Mirror the exact bin masking in fft_bandpass_filtering.m.
        outside = np.flatnonzero((bins < low) | (bins > high))
        spectrum = fft.copy()
        spectrum[outside] = spectrum[n - 1 - outside] = 0
        return np.fft.ifft(spectrum).real

    envelopes = np.stack([np.abs(signal.hilbert(band(f-1, f+1))) for f in freqs], axis=1)
    log_env = np.log10(np.maximum(envelopes, np.finfo(float).tiny))
    edges = np.floor(np.linspace(0, n-1, 5) + 0.5).astype(int)
    fits = [background(freqs, np.median(log_env[a:b+1], axis=0)) for a,b in zip(edges, edges[1:])]
    valid = [fit for fit in fits if fit is not None]
    if not valid:
        raise ValueError("No quarter of this recording provided enough bins for the background fit.")
    residual = log_env - min(valid, key=lambda fit: fit[0])
    threshold = 2 * residual.std(ddof=1)
    boxes = []
    for f in freqs:
        for a,b in runs(residual[:, f-1] > threshold):
            if f * (b-a) / fs < 2:
                continue
            hz = confirmed_frequency(raw[a:b+1], fs, f)
            if hz is not None:
                peak = residual[a:b+1, max(f-2,0):min(f+1,len(freqs))].max()
                boxes.append(Box(a,b,f-1,f+1,hz,peak))

    taper = signal.windows.tukey(n, 0.25)

    @lru_cache(maxsize=32)
    def tapered(low, high):
        return taper * band(low, high)

    i = 0
    while i < len(boxes):
        first = boxes[i]
        for j in range(i+1, len(boxes)):
            other = boxes[j]
            a, b = max(first.start, other.start), min(first.stop, other.stop)
            if b-a < 2 or abs(first.center-other.center) > 2:
                continue
            r, p = pearsonr(tapered(first.low,first.high)[a:b+1], tapered(other.low,other.high)[a:b+1])
            if r > 0 and p < 0.05:
                strongest = max((first,other), key=lambda box: box.peak)
                boxes[i] = Box(min(first.start,other.start), max(first.stop,other.stop),
                               min(first.low,other.low), max(first.high,other.high),
                               strongest.center, strongest.peak)
                boxes.pop(j)
                i = 0
                break
        else:
            i += 1
    return residual, boxes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parents[1]/"data"/"data_for_figure8.mat")
    parser.add_argument("--out", type=Path, default=Path("cho-demo-output"))
    args = parser.parse_args()
    raw = loadmat(args.input, variable_names=["raw_data"])["raw_data"].squeeze()
    fs = 500
    residual, boxes = detect(raw,fs)
    args.out.mkdir(exist_ok=True, parents=True)
    with (args.out/"events.csv").open("w",newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["start","stop","low","high","center","peak"])
        writer.writeheader()
        writer.writerows(asdict(box) for box in boxes)
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle
    fig = Figure(figsize=(10,5.5),layout="constrained")
    trace, tf = fig.subplots(2,1,sharex=True)
    t = np.arange(len(raw))/fs
    trace.plot(t,raw,color="#243e54",linewidth=0.8)
    trace.set(ylabel="Signal amplitude",title="Hippocampal recording")
    tf.imshow(residual.T, origin="lower", aspect="auto", extent=[-0.5/fs,t[-1]+0.5/fs,0.5,40.5],
              vmin=0,vmax=3*residual.std(),cmap="viridis")
    for box in boxes:
        tf.add_patch(Rectangle((box.start/fs,box.low),(box.stop-box.start)/fs,box.high-box.low,
                               fill=False,edgecolor="white",linewidth=1.2))
        tf.annotate(f"{box.center} Hz", ((box.start+box.stop)/(2*fs),box.high),
                    xytext=(0,9),textcoords="offset points",ha="center",va="bottom",
                    color="white",fontsize=10,weight="bold",
                    bbox=dict(boxstyle="round,pad=0.2",fc="#152438",alpha=.9,ec="none"))
    tf.set(xlabel="Time (s)",ylabel="Frequency (Hz)",title=f"Background-corrected envelope · {len(boxes)} detected episodes")
    fig.savefig(args.out/"detections.png",dpi=160)
    print(f"{len(raw)} samples at {fs} Hz; {len(boxes)} episodes. Output: {args.out}")


if __name__ == "__main__":
    main()
