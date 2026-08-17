"""Generate the EEG<->MNE numeric-parity reference for PhysioAnalysis (VAL-04).

Replaces the manuscript's documentation-level "MNE-equivalent" correspondence
with a bundled, runnable numeric reference captured from MNE-Python:
  * Welch PSD band power (MNE mne.time_frequency.psd_array_welch) for the
    PhysioAnalysis default EEG bands -> the manuscript's "alpha power" metric.
  * PLV (alpha) and wPLI (theta) from mne_connectivity.spectral_connectivity_epochs.

The SAME fixed multi-channel signal (seeded, with graded alpha coupling and a
lagged theta source) drives both MNE and PhysioAnalysis, so the R side runs
bandPower()/plv()/wPLI() on identical input. Writes /tmp/val04/eeg_mne_reference.json;
data-raw/eeg_mne_reference.R builds tests/testthat/fixtures/eeg-mne-reference.rds.

Note: PhysioAnalysis's plv()/wPLI() are continuous Hilbert estimators while
mne_connectivity uses an epoch-spectral estimator, so PLV/wPLI agree on coupling
STRUCTURE (~0.95 correlation) but are not identical; PSD band power matches to
< 0.5%. Run:  micromamba run -n mne python data-raw/eeg_mne_reference.py
"""
import json, os
import numpy as np
import mne
from mne.time_frequency import psd_array_welch
from mne_connectivity import spectral_connectivity_epochs
from scipy.signal import butter, filtfilt
mne.set_log_level("ERROR")

OUT = os.environ.get("OUT", "/tmp/val04")
os.makedirs(OUT, exist_ok=True)
sf, T, nch = 256.0, 20.0, 8
n = int(sf * T)
rng = np.random.default_rng(42)

# PhysioAnalysis bandPower() default bands
BANDS = {"delta": (0.5, 4), "theta": (4, 8), "alpha": (8, 13),
         "beta": (13, 30), "gamma": (30, 100)}

def pink(n, amp):
    x = np.fft.irfft(np.fft.rfft(rng.standard_normal(n)) /
                     np.sqrt(np.maximum(np.fft.rfftfreq(n, 1 / sf), 1 / n)), n)
    return x / np.std(x) * amp

def narrow(n, lo, hi, amp):
    b, a = butter(4, [lo / (sf / 2), hi / (sf / 2)], btype="band")
    y = filtfilt(b, a, rng.standard_normal(n))
    return y / np.std(y) * amp

data = np.vstack([pink(n, a) for a in np.linspace(3e-6, 6e-6, nch)])
common_a = narrow(n, 8, 12, 1.0)
w = np.array([0.9, 0.85, 0.75, 0.5, 0.2, 0.1, 0.1, 0.05])
for i, amp in enumerate(np.linspace(4e-6, 11e-6, nch)):
    ai = w[i] * common_a + (1 - w[i]) * narrow(n, 8, 12, 1.0)
    data[i] += amp * ai / np.std(ai)
common_th = narrow(n, 4, 8, 1.0)
for i, lag in zip(range(4), [0, 3, 6, 9]):
    data[i] += 5e-6 * np.roll(common_th, lag)

psd, freqs = psd_array_welch(data, sf, fmin=0, fmax=sf / 2, n_per_seg=256,
                             n_overlap=128, window="hann", verbose=False)
band_powers = {nm: psd[:, (freqs >= lo) & (freqs <= hi)].sum(axis=1).tolist()
               for nm, (lo, hi) in BANDS.items()}

ep_len = 512
n_ep = n // ep_len
epd = data[:, :n_ep * ep_len].reshape(nch, n_ep, ep_len).transpose(1, 0, 2)
def conn(method, lo, hi):
    c = spectral_connectivity_epochs(epd, method=method, sfreq=sf, fmin=lo, fmax=hi,
                                     faverage=True, verbose=False)
    m = c.get_data(output="dense")[:, :, 0]
    return (m + m.T).tolist()

json.dump(dict(sf=sf, nch=nch, n=n, bands=list(BANDS),
               data=np.round(data, 12).tolist(),
               band_powers=band_powers,
               plv=conn("plv", 8, 12), wpli=conn("wpli", 4, 8),
               mne=mne.__version__),
          open(os.path.join(OUT, "eeg_mne_reference.json"), "w"))
print("wrote", os.path.join(OUT, "eeg_mne_reference.json"), "| mne", mne.__version__)
