# PyAMDF.py, AMDF-demo.py, RFFT-demo.py
Python implementation of Average Magnitude Difference Function for speech F0 estimation. Uses builtins sys, re and numpy, scipy, scipy, matplotlib. PyAMDF.py provides a detailed set of plots illustrating various transformations of the speech signal. AMDF-demo.py provides a set of plots illustrating the moving window and subtraction operation. RFFT-demo.py is a demo of estimating the real component of the FFT, for the Long-Term spectrum (LTS) and for the Amplitude Envelope Modulation Spectrum (AEMS). The PNG file longtermspectrum-aeksjiayan-5s.png shows the kind of output provided. AEM is extracted with a peak maximum technique over a moving window, rather than Hilbert Transform.

The fundamental frequency (aka 'F0', 'pitch') of the voice is generated by air flowing through the vocal folds in the larynx. It is interrupted by pauses and by voiceless consonants and is absent in whispering. Where the fundamental frequency is not interrupted, it is a carrier signal, modulated by changing the air pressure during expiration and by changing the rigidity of the vocal folds. This frequency modulation (FM) is entirely analogous to the FM of radio broadcasts, except that the carrier signal of speech is between 50 and 500 Hz, while the carrier signal of FM radio is around 100 MHz, and that the frequency of the modulation in speech is around 3 frequency changes per second, while the modulation of FM radio is the whole of speech, music or other sounds.

The methods of estimating speech F0 fall into two main types, which can also be combined: time domain methods, in which similarities in amplitude patterns between intervals in the waveform are analysed, and frequency domain methods, in which frequency patterns in intervals in the waveform are analysed.

The PyAMDF method is a time-domain method, and uses the Average Magnitude Difference Function, which checks the differences (better: similarity) intervals in the waveform and following intervals of the same length. The interval between the beginning of the first interval and the beginning of the next most similar interval is the period estimation, the inverse of the fundamental frequency estimation.

The code is mainly in traditional functional style, so that those familiar with other programming languages should not have too much difficulty in following it. The structure is:

	- Import of supporting library modules
	
	- Definition of parameters
	
	- Inport of a wavefile and calculation of relevant dependent parameters
	
	- Preprocessing of the signal by centre-clipping
	
	- F0 estimation by AMDF
	
	- Postprocessing of the signal by
	
		- Truncating F0 estimates outside a defined voice range
		
		- Smoothing the F0 pattern
		
	- Display of the waveform, a spectrogram and the F0 estimation

The PyAMDF implementation does not include common features of well-known standard F0 estimators such as Praat, RAPT (aka 'get_f0' or 'esps'), Reaper, SWIPE, YAAPT or YIN (descriptions of these can be found on the internet):

	- separate voicing detection
	
	- statistical estimation of parameter costs
	
	- octave and suboctave detection
	
	- continuity detection
	
Consequently, PyAMDF is not as robust or versatile as the F0 estimators named above, but is neverthess useful with careful selection of parameter values.

The implementation is standalone code, intended for use on the command line interface (CLI). It was developed and tested with Ubuntu Linux 16.10, but is intended to be platform-independent, but requires installation of the standard Python library modules NumPy, MatPlotLib and SciPy.

Typical CLI usage of PyAMDFis:

	PyAMDF.py somebodysvoice.wav fairlyhigh single

CLI parameter shortcuts for different voice type parameters are provided:

	'high', 'fairlyhigh', 'mid', 'fairlylow', 'low'

Other parameters are 'hard-wired' in the code, and can be changed as desired.

Three outputs are provided:

	- a screen figure display of waveform, spectrogram and F0 estimation,
	
	- a PNG file with the figure,
	
	- a CSV file of the CSV estimation for input into a spreadsheet.

The code is open source free software. The only condition on using this code is that it may not be patented or copyrighted. If the code is used in publications it may be cited as:

Gibbon, Dafydd. 2018. PyAMDF: a do-it-yourself pitch estimator in Python. Bielefeld University.

https://github.com/dafyddg/PyAMDF/
