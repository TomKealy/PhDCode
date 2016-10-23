M = 300;
K = 100;

sigma_squared = [1, 10, 100, 1000];
noise = normrnd(0, 1, [1, M]);
runs = 1;% size(sigma_squared,2);
max_not_in_positions = 0;

positions = randi(M, [1, 5]);
weights = normrnd(100, sigma_squared(4), [1, 5]);
     
x = zeros(1,M);
x(positions) = weights;

xe = x + noise;

A = normrnd(0, 1/M, [K, M]);
 
[U, S, V] = svd(A);

y = A*real(fft(x'));

xbar = V'*real(fft(x'));
ybar = U'*y;

n = 1:M;
