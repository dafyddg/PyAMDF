#!/usr/bin/python
# RFFT-demo.py

import sys
import numpy as np
import scipy.io.wavfile as wave
import matplotlib.pyplot as plt

# Frequency scale is straightforward:
# simply set xlim to desired values in Hz

# Define parameters
wavfile = sys.argv[1]
wavfilebase = wavfile.split('.')[0]
fmin = int(sys.argv[2]); fmax = int(sys.argv[3])
aemsfmin = int(sys.argv[4]); aemsfmax = int(sys.argv[5])

# Define envelope maxpeak window
envwin = 5
# also variables for aemspower and aemds power!!!

# Input WAV file and calculate params
sample_rate,signal = wave.read(wavfile)
downsamp = 4
sample_rate = sample_rate/downsamp
signal = signal[::downsamp]

sample_period = 1.0/sample_rate
siglen = len(signal)
sigdur = 1.0*siglen/sample_rate

# Calculate long term FFT
fouriercoeff = abs(np.fft.rfft(signal))
fouriercoeff = 1.0*fouriercoeff / np.max(fouriercoeff)
frequencies = np.fft.rfftfreq(signal.size, d=sample_period)

# Calculate amplitude modulation envelope
abssig = np.abs(signal); abssiglen = len(abssig)
envrange = range(abssiglen-envwin)
envelope = [ np.max(abssig[i:i+envwin]) for i in envrange]+[np.median(abssig)]*envwin
envelope = np.array(envelope)

# Calculate AEMS
aemsfourcoeff = abs(np.fft.rfft(envelope))
aemsfourcoeff = 1.0*aemsfourcoeff / np.max(aemsfourcoeff)
aemsfourfreq = np.fft.rfftfreq(envelope.size,d=sample_period)
aemsfourcoeff[0] = np.median(aemsfourcoeff)

# CALCULATE AEMDS
aemsfourcoeffdiff = np.append([0],np.abs(np.diff(aemsfourcoeff)))
aemsfourcoeffdiff[0] = 0

#=====================================================================
# Figure
_,(sig,lts,sigenv,aem,aems,aemds) = plt.subplots(6,1,figsize=(14,12))

# Waveform
sig.set_title("Waveform")
sig.set_xlabel( "Time (s), samprate=%d,siglen=%d,sigdur=%.3f"%(sample_rate,siglen,sigdur) )
sig.set_xlim(0,sigdur)
sig.set_yticks([])
x = np.linspace(0,siglen,siglen)/sample_rate
# x = range(len(signal))
sig.plot(x,signal,color='green')

# Long term spectrum
lts.set_title("Long term spectrum (frequency x fourier coefficients)")
lts.set_xlabel("Frequency (Hz)")
lts.set_xlim(fmin,fmax)	# Max freq = len(frequencies) / 2
lts.set_ylim(0,1)
lts.plot(frequencies,fouriercoeff)

# Waveform
sigenv.set_title("Waveform, Rectified Waveform, Amplitude Envelope Modulation spectrum")
sigenv.set_xlabel( "Time (s)")
sigenv.set_xlim(0,sigdur)
sigenv.set_yticks([])
x = np.linspace(0,siglen,siglen)/sample_rate
# x = range(len(signal))
sigenv.plot(x,signal,color='green')
sigenv.plot(x,np.abs(signal),color='lightgreen')
sigenv.plot(x,envelope,linewidth=3,color='r')

# Amplitude envelope
aem.set_title("Amplitude Modulation Envelope")
aem.set_xlabel("Time (s)")
aem.set_yticks([])
x = np.linspace(0,siglen,siglen)/sample_rate
aem.plot(x,envelope,color='r')

# Amplitude envelope magnitude spectrum
aems.set_title("Amplitude Envelope Modulation Spectrum")
aems.set_xlabel("Frequency (Hz)")
aems.set_xlim(aemsfmin,aemsfmax)
aems.stem(aemsfourfreq,aemsfourcoeff/np.max(aemsfourcoeff))

# Amplitude envelope magnitude spectrum
aemds.set_title("Amplitude Envelope Modulation Difference Spectrum")
aemds.set_xlabel("Frequency (Hz)")
aemds.set_xlim(aemsfmin,aemsfmax)
aemds.stem(aemsfourfreq, aemsfourcoeffdiff/np.max(aemsfourcoeffdiff), color='purple')

plt.tight_layout(pad=1,w_pad=0.5,h_pad=0.5)
plt.savefig("longtermspectrum-aems%s.png"%wavfilebase)
plt.show()

