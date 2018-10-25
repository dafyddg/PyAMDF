#!/usr/bin/python
# spectrogram4laura.py
# D. Gibbon 2018-09-03

# Usage: spectrogram4laura.py <yourfile.wav>

try:
	import sys
	import matplotlib.pyplot as plt
	import scipy.io.wavfile as wave
	from scipy.signal import spectrogram
except:
	print "The library modules sys, matplotlib, scipy must be installed."
	exit()
try:
	sampfreq,signal = wave.read(sys.argv[1])
	plt.specgram(signal, NFFT=1000, Fs=sampfreq)
	plt.ylim(0,8000)
	plt.show()
except:
	print "Usage: spectrogram yourwavefile.wav (must be mono)"
	exit()
