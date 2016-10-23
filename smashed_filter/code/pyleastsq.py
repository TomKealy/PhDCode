import numpy as np
import matplotlib.pyplot as plt
import scipy

# set up true function and "measured" data
x = np.linspace(0, 18e-2, 300);
A, k, theta = 10, 1.0 / 3e-2, np.pi / 6;
y_true = A * np.sin(2 * np.pi * k * x + theta);
y_meas = y_true + 2*np.random.randn(x.size);
plt.plot(x, y_meas);
plt.plot(x, y_true);
plt.show()

# residual function, e_i
def residuals(p, y, x):
    A, k, theta = p;
    err = y - A * np.sin(2 * np.pi * k * x + theta);
    return err;

def peval(x, p):
    return p[0] * np.sin(2 * np.pi * p[1] * x + p[2]);

# starting values of A, k and theta
p0=[12.0, 35.0, np.pi]
print(np.array(p0));

# do least squares fitting
from scipy.optimize import leastsq

plsq = leastsq(residuals, p0, args=(y_meas, x));
print(plsq[0]); print(np.array([A, k, theta]));

plt.plot(x, peval(x, plsq[0]))
plt.plot(x, y_meas,'ro')
plt.plot(x, y_true);
plt.title('Least-squares fit to noisy data');
plt.legend(['Fit', 'Noisy', 'True']);
