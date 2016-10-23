function y = guassianpulse(rate, decay, bandwidth)

t = 0 : 1/rate : 10e-3;
d = [0 : 1/1e3 : 10e-3 ; decay.^(0:10)]';
y = pulstran(t,d,'gauspuls',10e3,0.7);

plot(t,y)
xlabel 'Time (s)', ylabel 'Periodic Gaussian pulse'

end