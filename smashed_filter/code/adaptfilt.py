import numpy as np
import matplotlib.pyplot as plt
import adaptfilt as adf

u = np.load('speech.npy')

coeffs = np.concatenate(([0.8], np.zeros(8), [-0.7], np.zeros(9), [0.5], np.zeros(11), [-0.3], np.zeros(3), [0.1], np.zeros(20), [-0.05]  ))

d = np.convolve(u, coeffs)

v = np.random.randn(len(d))*np.sqrt(5000)

d += v

M = 100

step = 0.1

y, e, w, = adf.nlms(u, d, M, step, returnCoeffs=True)

mswe = adf.mswe(w, coeffs)

# Plot speech signals
plt.figure()
plt.title("Speech signals")
plt.plot(u, label="Emily's speech signal, u(n)")
plt.plot(d, label="Speech signal from John, d(n)")
plt.grid()
plt.legend()
plt.xlabel('Samples')

