ECG Signal Filtering (FIR vs IIR)
Removes the three most common noise sources from an ECG and compares an IIR (Butterworth + notch) pipeline with a FIR (Hamming-window) one.

Noise	Cause	Fix
Baseline wander (~0.3 Hz)	Breathing, movement	High-pass at 0.5 Hz
Powerline interference (50 Hz)	Mains hum	Notch / low-pass
Muscle (EMG) noise	Muscle activity	Low-pass at 40 Hz
A synthetic ECG (P, Q, R, S, T Gaussian waves, 72 bpm, fs = 360 Hz) is used so the clean reference is known and SNR can be measured exactly. The same code works on real recordings (e.g. MIT-BIH from PhysioNet); load the samples into noisy instead.

Run
pip install -r requirements.txt
python src/ecg_filter.py
Results
Signal	Output SNR
Noisy	-8.20 dB
IIR (Butterworth 0.5-40 Hz + 50 Hz notch)	11.87 dB
FIR (1001-tap Hamming 0.5-40 Hz)	11.52 dB
Time domain Spectrum Filter responses

Takeaways
The IIR filter reaches the same quality with an order-4 filter instead of 1001 taps, which is far cheaper computationally.
The FIR filter has exactly linear phase, which preserves the shape of the QRS complex. Both are applied with filtfilt (zero-phase), so neither distorts timing here.
A real-time version would need causal filtering (lfilter/sosfilt); the FIR would then add a fixed delay, while the IIR would introduce phase distortion.
Ideas to extend
Load a real MIT-BIH record with the wfdb package.
Add R-peak detection (Pan-Tompkins) and compute heart rate.
Implement the filters in C for a microcontroller.
‎01-ecg-signal-filtering/requirements.txt‎
+3
Lines changed: 3 additions & 0 deletions
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,3 @@
numpy
scipy
matplotlib
‎01-ecg-signal-filtering/results/filter_response.png‎
60.3 KB

‎01-ecg-signal-filtering/results/spectrum.png‎
81.7 KB

‎01-ecg-signal-filtering/results/time_domain.png‎
275 KB

‎01-ecg-signal-filtering/src/ecg_filter.py‎
+132
Lines changed: 132 additions & 0 deletions
Original file line number	Diff line number	Diff line change
@@ -0,0 +1,132 @@
"""
ECG Signal Filtering using FIR and IIR filters.
Generates a synthetic ECG, corrupts it with the three classic noise sources
(baseline wander, 50 Hz powerline interference, high-frequency EMG noise),
then cleans it with:
* IIR: Butterworth band-pass + notch filter (zero-phase via filtfilt)
* FIR: Hamming-window band-pass (zero-phase via filtfilt)
and compares them using output SNR.
Run:  python ecg_filter.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal
FS = 360          # sampling rate (Hz), same as MIT-BIH database
DURATION = 10     # seconds
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)
def synthetic_ecg(fs=FS, duration=DURATION, hr=72):
"""Build an ECG-like waveform from Gaussian P, Q, R, S, T waves."""
t = np.arange(0, duration, 1 / fs)
beat = 60.0 / hr
ecg = np.zeros_like(t)
# (centre offset as fraction of beat, width in s, amplitude)
waves = [(0.20, 0.025, 0.15),   # P
(0.37, 0.010, -0.15),  # Q
(0.40, 0.012, 1.00),   # R
(0.43, 0.010, -0.25),  # S
(0.65, 0.045, 0.30)]   # T
for k in range(int(duration / beat) + 1):
for c, w, a in waves:
ecg += a * np.exp(-((t - (k * beat + c * beat)) ** 2) / (2 * w ** 2))
return t, ecg
def add_noise(t, ecg, seed=1):
rng = np.random.default_rng(seed)
wander = 0.5 * np.sin(2 * np.pi * 0.3 * t)          # breathing / movement
powerline = 0.25 * np.sin(2 * np.pi * 50 * t)       # mains hum
emg = 0.08 * rng.standard_normal(len(t))            # muscle noise
return ecg + wander + powerline + emg
def iir_filter(x, fs=FS):
sos = signal.butter(4, [0.5, 40], btype="bandpass", fs=fs, output="sos")
y = signal.sosfiltfilt(sos, x)
b, a = signal.iirnotch(50, Q=30, fs=fs)
return signal.filtfilt(b, a, y)
def fir_filter(x, fs=FS):
taps = signal.firwin(1001, [0.5, 40], pass_zero=False, window="hamming", fs=fs)
return signal.filtfilt(taps, [1.0], x)
def snr_db(clean, test):
noise = test - clean
return 10 * np.log10(np.sum(clean  2) / np.sum(noise  2))
def main():
t, clean = synthetic_ecg()
noisy = add_noise(t, clean)
iir = iir_filter(noisy)
fir = fir_filter(noisy)
m = slice(FS // 2, -FS // 2)   # ignore edge transients
def prep(x): return x[m] - np.mean(x[m])
results = {
"Noisy": snr_db(prep(clean), prep(noisy)),
"IIR (Butterworth + notch)": snr_db(prep(clean), prep(iir)),
"FIR (Hamming)": snr_db(prep(clean), prep(fir)),
}
print("Output SNR (dB):")
for k, v in results.items():
print(f"  {k:28s} {v:6.2f}")
fig, ax = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
for a, y, title in zip(ax, [clean, noisy, iir, fir],
["Clean (reference)", "Noisy ECG", "IIR filtered", "FIR filtered"]):
a.plot(t, y, lw=0.9)
a.set_title(title, loc="left", fontsize=10)
a.set_ylabel("mV")
a.grid(alpha=0.3)
ax[-1].set_xlabel("Time (s)")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "time_domain.png"), dpi=130)
plt.close()
fig, ax = plt.subplots(figsize=(10, 4.5))
for y, lab in [(noisy, "Noisy"), (iir, "IIR"), (fir, "FIR")]:
f, p = signal.welch(y, FS, nperseg=1024)
ax.semilogy(f, p, label=lab)
ax.set_xlim(0, 100)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("PSD")
ax.set_title("Power spectral density before / after filtering")
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "spectrum.png"), dpi=130)
plt.close()
sos = signal.butter(4, [0.5, 40], btype="bandpass", fs=FS, output="sos")
w, h_iir = signal.sosfreqz(sos, worN=4096, fs=FS)
taps = signal.firwin(1001, [0.5, 40], pass_zero=False, window="hamming", fs=FS)
w2, h_fir = signal.freqz(taps, worN=4096, fs=FS)
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.plot(w, 20 * np.log10(np.abs(h_iir) + 1e-12), label="IIR Butterworth (order 4)")
ax.plot(w2, 20 * np.log10(np.abs(h_fir) + 1e-12), label="FIR Hamming (1001 taps)")
ax.set_ylim(-80, 5)
ax.set_xlim(0, 100)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Magnitude (dB)")
ax.set_title("Filter magnitude responses")
ax.grid(alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT, "filter_response.png"), dpi=130)
plt.close()
if __name__ == "__main__":
main()
