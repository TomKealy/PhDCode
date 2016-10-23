# Import the plotting library
from matplotlib import pyplot as plt
import scipy.signal
import numpy as np

# Define the setup
fCarrier = 10;
fAudio = 1;
fs = 1000;
timeEnd = 1;
time = np.linspace(0,2,fs*timeEnd);

# Create the signals
carrier = np.sin(2*np.pi*fCarrier*time);
audio = np.sin(2*np.pi*fAudio*time);
audioInt = -np.cos(2*np.pi*fAudio*time);
freqMod = np.sin(2*np.pi*fCarrier*time + 2*np.pi*1*audioInt);

plt.plot(freqMod)
plt.show()

# Downconvert
analyticSignal = scipy.signal.hilbert(freqMod); # wikipedia analytic signal
baseband = analyticSignal * np.exp(-2*np.pi*fCarrier*time*1j); # complex mixing
audioDemod = np.angle( baseband[1::1] * np.conjugate(baseband[0:-1:1]) ); # fm demo

plt.plot(analyticSignal)
plt.show()
