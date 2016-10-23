% Avetis Ioannisyan
% avetis@60ateight.com
% Last Updated: 11/30/05
% LMS Channel Adaptation

% reset randomizers
randn('state',sum(100*clock)) ;
rand('state',sum(100*clock)) ;

numPoints = 5000;
numTaps = 10;		% channel order
Mu = 0.01;          % iteration step size

% input is guassian
x = randn(numPoints,1) + 1j*randn(numPoints,1); 

% choose channel to be random uniform
h = rand(numTaps, 1);% + i*rand(numTaps, 1);

%h = [1 0 0 0 1];   %testing only
h = h/max(h);   %normalize channel

% convolve channel with the input
d = filter(h, 1, x);

% initialize variables
w = [];
y = [];
in = []; 
e = []; % error, final result to be computed

w = zeros(numTaps+1,1) + 1j*zeros(numTaps+1,1);

% LMS Adaptation
for n  = numTaps+1 : numPoints
    
    % select part of training input
    in = x(n : -1 : n-numTaps) ;

    y(n) = w'*in;        

    % compute error
    e(n) = d(n)-y(n);

    % update taps
    w = w + Mu*( real(e(n)*conj(in)) - 1j*imag(e(n)*conj(in)) );

end

[erls, w] = rls(0.9, numTaps, x, d, 1);

% Plot results
figure(10);
semilogy(1:numPoints, abs(e), 'b', 1:numPoints, abs(erls), 'r');
title(['LMS Adaptation Learning Curve Using Mu = ', num2str(Mu)]);
xlabel('Iteration Number');
ylabel('Output Estimation Error in dB');