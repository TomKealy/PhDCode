import numpy as np
from sklearn.gaussian_process import GaussianProcess
from matplotlib import pyplot as pl

np.random.seed(1)


def f(x):
    """The function to predict"""
    return x * np.sin(x)

X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T

y = f(X).ravel()

x = np.atleast_2d(np.linspace(0, 10, 1000)).T

gp = GaussianProcess(corr='cubic', theta0=1e-2, thetaL=1e-4, thetaU=1e-1, random_start=100)

gp.fit(X, y)

y_pred, MSE = gp.predict(x, eval_MSE=True)
sigma = np.sqrt(MSE)

fig = pl.figure()
pl.plot(x, f(x), 'r:', label=u"$f(x)=x\,\sin(x)")
pl.plot(X, y, 'r.', markersize=10, label=u"Observations")
pl.plot(x, y_pred, 'b-', label=u"Prediction")
pl.fill(np.concatenate([x, x[::-1]]),
        np.concatenate([y_pred - 1.9600 * sigma,
                       (y_pred + 1.9600 * sigma)[::-1]]),
        alpha=.5, fc='b', ec='None', label='95% confidence interval')

pl.xlabel("$x$")
pl.ylabel("$f(x)$")
pl.ylim(-10, 20)
pl.legend(loc="upper left")
pl.show()
