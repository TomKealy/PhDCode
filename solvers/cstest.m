clear all;

lena = imread('lena.jpg');

X = zeros(1,512);

positions = randi(512,[1,30]);

X(positions) = 1;

dctX = dct(X);

dctlena = dct(lena);

B = binornd(1,0.5,150,length(X));

y = B*X';

z = B*dctlena;

z0 = B'*z

x0 = B'*y;

xhat = l1eq_pd(x0, B, [], y, 1e-3);

zhat = l1eq_pd(z0, B, [], z, 1e-3);