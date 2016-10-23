import numpy
import pylab

# intial parameters
n_iter = 500
sz = (n_iter,) # size of array
x = numpy.array( [-0.37727 + j * j * 0.00001 for j in xrange(n_iter)]) # truth value (typo in example at top of p. 13 calls this z)
z = numpy.array([numpy.random.normal() * 2.0 - 1.0 + actual_val for actual_val in x])

#numpy.random.normal(x,0.1,size=sz) # observations (normal about x, sigma=0.1)

Q = 1e-5 # process variance

# allocate space for arrays
xhat=numpy.zeros(sz)      # a posteri estimate of x
P=numpy.zeros(sz)         # a posteri error estimate
xhatminus=numpy.zeros(sz) # a priori estimate of x
Pminus=numpy.zeros(sz)    # a priori error estimate
K=numpy.zeros(sz)         # gain or blending factor
posteri_estimate_for_graphing = []

R = 0.1**2 # estimate of measurement variance, change to see effect

# intial guesses
xhat[0] = 0.0
P[0] = 1.0

# allocate space for arrays
posteri_estimate_for_graphing = []

# intial guesses
posteri_estimate = 0.0
posteri_error_estimate = 1.0

for iteration in range(1, n_iter):
    # time update
    priori_estimate = posteri_estimate
    priori_error_estimate = posteri_error_estimate + Q

    # measurement update
    blending_factor = priori_error_estimate / (priori_error_estimate + R)
    posteri_estimate = priori_estimate + blending_factor * (z[iteration] - priori_estimate)
    posteri_error_estimate = (1 - blending_factor) * priori_error_estimate
    posteri_estimate_for_graphing.append(posteri_estimate)

for k in range(1,n_iter):
    # time update
    xhatminus[k] = xhat[k-1]
    Pminus[k] = P[k-1]+Q

    # measurement update
    K[k] = Pminus[k]/( Pminus[k]+R )
    xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
    P[k] = (1-K[k])*Pminus[k]

pylab.figure()
pylab.plot(z,'r',label='noisy measurements')
pylab.plot(xhat,'b-',label='a posteri estimate')
pylab.plot(x,color='g',label='truth value')
pylab.legend()
pylab.xlabel('Iteration')
pylab.ylabel('Voltage')

pylab.figure()
valid_iter = range(1,n_iter) # Pminus not valid at step 0
pylab.plot(valid_iter,Pminus[valid_iter],label='a priori error estimate')
pylab.xlabel('Iteration')
pylab.ylabel('$(Voltage)^2$')
pylab.setp(pylab.gca(),'ylim',[0,.01])

pylab.figure()
pylab.plot(z, color='r', label='noisy measurements')
pylab.plot(numpy.array(posteri_estimate_for_graphing), 'b', label='a posteri estimate')
pylab.plot(x, color='g', label='truth value')
pylab.legend()
pylab.xlabel('Iteration')
pylab.ylabel('Voltage')
pylab.show()
