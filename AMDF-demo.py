#!/usr/bin/python
# -*- coding: utf-8 -*-'

# AMDF-demo.py

# D. Gibbon
# 2018-10-10
# Pitch Estimator based on Average Magnitude Difference Function

import sys, re
import numpy as np
from scipy.signal import medfilt, hilbert
import scipy.io.wavfile as wav
import matplotlib.pyplot as plt

#==============================================================================
# Parameters

try:
	filename = sys.argv[1]
	voice = sys.argv[2]
	figtype = sys.argv[3]
	filebase = re.sub('.wav','',filename)
	print filename,filebase

except:
	print 'Usage: PyADMF.py wavfilename high|fairlyhigh|mid|fairlylow|low single|multiple'; exit()

try:
	if figtype == 'single':
		singlefigure = True
	else:
		singlefigure = False
		
	framerate = 0.01
	
	centreclip = 10
	winoffsetdivisor = 20.0
	medianwindow = 3

	envelopeflag = False
	
	spectmin = 0
	spectmax = 5000
	spectwinfactor = 20

	if voice == 'high':
		f0min = 160
		f0max = 450

	elif voice == 'fairlyhigh':
		f0min = 140
		f0max = 400

	elif voice == 'mid':
		f0min = 120
		f0max = 350

	elif voice == 'fairlylow':
		f0min = 100
		f0max = 300

	elif voice == 'low':
		f0min = 50
		f0max = 250

	else:
		print 'Unknown voice type.'; exit()

except:
	print 'Parameter error.'; exit()

#==============================================================================
#==============================================================================
# WAV file input

try:
	fs, signal = wav.read(filename)

	sigstart = 6400
	sigend = sigstart + 4/framerate		# for a good chunk of jiayan.wav

#	plt.plot(signal[sigstart:sigend]) ; plt.show() ; exit()	# for chunk test

	signallen = len(signal)
	signalduration = int(round(1.0*signallen/fs))

	sampframeratio = fs * framerate

	framelen = np.int(round(fs * framerate))
	framecount = np.int(round(1.0 * signallen / framelen))
	diffoffset = np.int(round(framelen/winoffsetdivisor))

	newsignallen = framecount * framelen
	signal = signal[:newsignallen]

	spectwin = int(fs / spectwinfactor)

	print "Sampling rate: %s Hz"%fs
	print "len(signal): %s"%len(signal)

except:
	print "Error reading signal."; exit()

#==============================================================================
#==============================================================================
# Preprocessing (centre-clipping)

try:
	"""
	signal = (abs(signal) > centreclip).astype(np.int) * signal
	"""

except:
	print "Error preprocessing."; exit()

#==============================================================================
# F0 estimation

if True:
	irange = range(0, signallen-2*framelen, framelen)	# Note truncation.
# Allocate memory for f0list and AMDF list
	f0list = np.zeros(framecount)
	meandiffs = np.zeros(framelen).tolist()
# Move frame window through signal

	smallestlist = []
	count = 0
	for framestart in irange:
		framestop = framestart + framelen
		frame = signal[framestart:framestop]
# Calculate Average Magnitude Difference Function with moving window
		indx = 0
		for winstart in range(framestart,framestop):
			movingwin = signal[winstart:winstart+framelen]
			absdiffs = abs(frame - movingwin)
			meandiffs[indx] = np.mean(absdiffs)
			indx += 1

# Pick smallest absolute difference in frame
		smallest = np.min(meandiffs[diffoffset:])
		
		smallestlist += [smallest]
# Get position of the smallest absolute difference

		index = meandiffs.index(smallest)
# Divide the sampling rate by the number of samples in the interval
		f0 = fs / index	# That is: t = index/fs; f0 = 1/t
# Extend f0 list
		f0list[count] = f0
		count += 1
		
		teststart = sigstart
		
		if framestart == teststart:
			segmentsize = 3 * framelen

			copystart = teststart + index
			framecopy = signal[copystart:copystart+framelen]

			absdiffs = abs(frame-framecopy)
			absdiffs = absdiffs/float(max(absdiffs))

			copylen = len(framecopy)
			signal = signal/float(max(abs(signal)))
			framecopy = framecopy/float(max(abs(framecopy)))
		
			period = float(index)/fs
			frequency = fs/index

			sigsegment = np.array(signal[teststart:teststart+segmentsize])
			sigsegmentlen = float(len(sigsegment))
			x = np.linspace(0.0,sigsegmentlen,sigsegmentlen)
			xx = np.array(x)/fs
			
			fontsize = 20

#==============================================================================

			fig, (sp1,sp2,sp3,sp4) = plt.subplots(4, 1, figsize=(22,8))

#			sp1.set_title("Smallest Average Magnitude Difference: Position: %s samples, Period: %.4f s, Frequency: %d Hz."%(index,period,frequency),fontsize=fontsize)
			sp1.plot(xx,sigsegment,label="Signal",linewidth=2)
			sp1.axvline(x=0.00003, label="Frame limits",linewidth=2,color='r')
			for i in np.linspace(0,sigsegmentlen,framelen/40)/fs:
				sp1.axvline(x=i,linewidth=2,color='r')
			sp1.set_xlabel("Time (s), sampling rate = %d Hz, segment length = %.3f s, frame length = %3f s"%(fs,1.0*sigsegmentlen/fs,1.0*framelen/fs),fontsize=fontsize)
			sp1.legend()

			sp2.plot(signal[teststart:teststart+segmentsize],label="Signal")
			sp2.plot(range(framelen),signal[teststart:teststart+framelen],color='r')
			sp2.scatter(range(framelen),signal[teststart:teststart+framelen],color='r',s=10,label="Frame")
			sp2.legend()
			for i in np.linspace(0.5,sigsegmentlen,framelen/40):
				sp2.axvline(x=i,linewidth=2,color='r')
			sp2.set_xlabel("Samples, sampling rate = %d Hz, sample period = %7f s"%(fs, 1.0/fs),fontsize=fontsize)

			sp3.plot(signal[teststart:teststart+segmentsize],label="Signal")
			sp3.plot(range(index,index+framelen),signal[copystart:copystart+framelen],color='g',label="Copy")
			sp3.scatter(range(index,index+len(framecopy)),framecopy,marker="x",color='g',s=30)
			for i in np.linspace(0.5,sigsegmentlen,framelen/40):
				sp3.axvline(x=i,linewidth=2,color='r')
			sp3.axvline(x=index,linewidth=4,color='purple',label="First smallest average difference")
			sp3.legend()
			
			sp4.plot(range(index,index+framelen),absdiffs,label="AMDF",color='orange',linewidth=2)
			sp4.axvline(x=index,linewidth=4,color='purple',label="First smallest average difference")
			sp4.axhline(xmin=0.0,xmax=index/sigsegmentlen, y=0.5,linestyle=':', linewidth=4, color='purple', label="Period: %.6fs, Frequency: %d Hz"%(1.0*index/fs,fs/index))
			for i in np.linspace(0.5,sigsegmentlen,framelen/40):
				sp4.axvline(x=i,linewidth=2,color='r')
			sp4.legend()

			sp1.set_xlim(0,1.0*sigsegmentlen/fs); sp1.set_ylim(-1,1)
			sp2.set_xlim(0,segmentsize); sp2.set_ylim(-1,1)
			sp3.set_xlim(0,segmentsize); sp3.set_ylim(-1,1)
			sp4.set_xlim(0,segmentsize); sp4.set_ylim(0,1)

if False:
	print "F0 estimation error."; exit()

#==============================================================================
# Post-processing

try:
# Remove f0 values outside defined limits
#	f0list = (f0list > f0min).astype(int) * f0list
#	f0list = (f0list < f0max).astype(int) * f0list
# Smooth F0 contour
#	f0list = medfilt(f0list,medianwindow)

#	plt.scatter(range(len(signal)),signal)
#	plt.scatter(range(len(smallestlist)),smallestlist)
#	plt.scatter(range(len(f0list)),f0list)

	plt.tight_layout()
	plt.savefig("speech-amd-frames-%s.png"%filebase)
	plt.show()

except:
	print "Post-processing failed."; exit()

#==============================================================================
#==============================================================================

exit()

#==============================================================================
#==============================================================================
# Print to CSV file (data in rows)

try:
	csvfilename = 'speech-f0estimate-%s.csv'%filebase
	f0csv = filebase+','+','.join(map(str,f0list))
	
	header = 'speech-'+filebase+','+','.join(map(str,np.linspace(0,len(f0list),len(f0list))*framerate))
	csvfiletext = header+'\n'+f0csv+'\n'
	file = open(csvfilename,'w')
	file.write(csvfiletext)
	print 'Output in CSV file. Open with "soffice %s", choose comma separator, click on the left of row 1, then Insert, then Chart, then format the chart as a scatter plot.'%csvfilename

except:
	print "Text file output error."; exit()

#==============================================================================
#==============================================================================
# Figures

try:
	sp1figname = 'speech-waveform-%s.png'%filebase
	sp2figname = 'speech-amenvelope-%s.png'%filebase
	sp3figname = 'speech-f0track-%s.png'%filebase
	sp4figname = 'speech-FMandAMenvelopes-%s.png'%filebase
	sp5figname = 'speech-spectrogram-%s.png'%filebase
	spfullfigname = 'speech-full-%s.png'%filebase

	spfignames = csvfilename,sp1figname,sp4figname,sp5figname,spfullfigname
	print "Output file names:\n\t%s\n\t%s\n\t%s\n\t%s\n\t%s"%spfignames

	if singlefigure:
		fig, (sp1,sp2,sp3,sp4,sp5) = plt.subplots(5, 1, figsize=(12,9))
		fig.suptitle('AMDF F0 estimation, DG, 2010-10-10',fontsize=16)
except:
	""; exit()

#==============================================================================
# Plot signal

if True:
	if not singlefigure:
		fig, (sp1) = plt.subplots(1, 1, figsize=(16,4))

	x = np.linspace(0,len(signal),len(signal))/fs
	signal = medfilt(signal,31)
	signal = signal / float(np.max(abs(signal)))
	abssignal = abs(signal)
	sp1.scatter(x,signal,color='g',s=2)
	if envelopeflag:
		sp1.plot(x,abssignal,color='lightgreen')
		peakwin = 20
		peakarray = [ max(abssignal[i-peakwin:i]) for i in range(peakwin,len(abssignal)) ]
		peakarray = np.append(medfilt(peakarray,201),[0]*peakwin)
		peakarray = peakarray / max(peakarray)
		sp1.plot(x,peakarray,color='r',linewidth=2)
	sp1.set_title('Waveform')
	sp1.set_xlabel('Time (s)')
	signalduration = float(len(signal))/fs
	sp1.set_xlim(0,signalduration)

	if not singlefigure:
		plt.tight_layout(pad=1, w_pad=0.1, h_pad=0.1)
		plt.savefig(sp1figname)

#==============================================================================
# Plot AM envelope

if True:
	if not singlefigure:
		fig, (sp2) = plt.subplots(1, 1, figsize=(12,4))
	x = np.linspace(0,len(signal),len(signal))/fs
	signal = medfilt(signal,31)
	signal = signal / float(np.max(abs(signal)))
	abssignal = abs(signal)
	peakwin = 20
	amenvelope = [ max(abssignal[i-peakwin:i]) for i in range(peakwin,len(abssignal)) ]
	amenvelope = np.append(medfilt(amenvelope,201),[0]*peakwin)
	amenvelope = amenvelope / max(amenvelope)
	sp2.plot(x,amenvelope,color='r',linewidth=2)
	sp2.set_title('Amplitude Modulation (AM) envelope')
	sp2.set_xlabel('Time (s)')
#	sp2.set_xlim(0,signalduration)
	if not singlefigure:
		plt.tight_layout(pad=1, w_pad=0.1, h_pad=0.1)
		plt.savefig(sp2figname)

#==============================================================================
# F0 estimation
try:
	if not singlefigure:
		fig, (sp3) = plt.subplots(1, 1, figsize=(12,4))
	sp3.set_title('AMDF Frequency Modulation (FM) envelope (F0 estimation, pitch tracking)')
	sp3.set_xlabel('Time (s)')
	sp3.set_xlim(0,signalduration)
	sp3.set_ylim(f0min,f0max)
	y = f0list
	leny = len(y)
	x = np.linspace(0,leny,leny)*framerate
	sp3.scatter(x,y,s=5,color='b')
	if not singlefigure:
		plt.tight_layout(pad=1, w_pad=0.1, h_pad=0.1)
		plt.savefig(sp3figname)
except:
	print "Problem with F0 estimation."; exit()

#==============================================================================
# F0 estimation with Average Magnitude Difference Function in moving window
try:
	if not singlefigure:
		fig, (sp4) = plt.subplots(1, 1, figsize=(12,4))
	if envelopeflag:
		amstring = ' AM envelope & '
	else:
		amstring = ''
	sp4.set_title(amstring+'AMDF Frequency Modulation (FM) envelope (F0 estimation, pitch tracking)')
	sp4.set_xlabel('Time (s)')
	sp4.set_xlim(0,signalduration)
	sp4.set_ylim(f0min,f0max)
	y = f0list
	leny = len(y)
	x = np.linspace(0,leny,leny)*framerate
	if envelopeflag:
		thinenvelope = amenvelope[::int(sampframeratio)]
		y = 0.5*y/max(y)
		z = thinenvelope/max(thinenvelope)
		sp4.plot(x,z,linewidth=2,color='r')
		sp4.set_ylim(0.0,1.1)
		sp4.set_yticklabels([])
	sp4.scatter(x,y,s=5,color='b')
	if not singlefigure:
		plt.tight_layout(pad=1, w_pad=0.1, h_pad=0.1)
		plt.savefig(sp4figname)

#==============================================================================
# Plot spectrogram

	if not singlefigure:
		fig, (sp5) = plt.subplots(1, 1, figsize=(12,8))
	sp5.set_title('Spectrogram')
	sp5.set_xlabel('Time (s)')
	sp5.specgram(signal,NFFT=spectwin, Fs=fs)
	sp5.axis(ymin=spectmin, ymax=spectmax)
	sp5.grid(which='both',axis='both',linewidth="1",linestyle='--')
	sp5.set_xlim(0,signalduration)
	if not singlefigure:
		plt.tight_layout(pad=1, w_pad=0.1, h_pad=0.1)
		plt.savefig(sp5figname)

#==============================================================================
# Format and display figure
	if singlefigure:

		plt.tight_layout(pad=3, w_pad=0.1, h_pad=1)
		plt.savefig(spfullfigname)
		plt.show()

except:
	print "Graphics or CSV output error."; exit()
#===========================================================================EOF
