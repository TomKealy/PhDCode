% EXAMPLE5 Polyphase filter implementation.

% This file is part of the examples provided for the MIT course
% 18.327 / 1.130 "Wavelets, Filter Banks and Applications"
% (http://web.mit.edu/1.130).
	
% Author : Kevin Amaratunga, MIT
% Version: 1.0


% Determine the filters.
clear all;
[h0,h1,f0,f1] = orthfilt(dbwavf('db2'));
M = length(h0);

% Load a signal.
load noisdopp
x = noisdopp;
N = length(x);

% Change the signal into polyphase form.
xeven = dyaddown(x,1);  % Even part
xodd = dyaddown(x,0);   % Odd part
xodd = [0 xodd(1:N/2-1)];
X = [xeven; xodd];

% Construct the polyphase matrix.
H = zeros(2,2,M/2);
H(1,1,:) = dyaddown(h0,1);  % h0,even[n]
H(1,2,:) = dyaddown(h0,0);  % h0,odd[n]
H(2,1,:) = dyaddown(h1,1);  % h1,even[n]
H(2,2,:) = dyaddown(h1,0);  % h1,odd[n]

% Run the polyphase filter.
Y = polyfilt(H,X);

% Plot the results.
L = N/2;
n = 0:L-1;
clf
subplot(2,1,1)
plot(n,Y(1,1:L))
axis tight
xlabel('Sample number')
ylabel('Lowpass')
title('Output from polyphase filter')
subplot(2,1,2)
plot(n,Y(2,1:L))
axis tight
xlabel('Sample number')
ylabel('Highpass')
pause

% Compute the results using the direct approach.
y0 = dyaddown(conv(x,h0),1);
y1 = dyaddown(conv(x,h1),1);

% Now compare the results.
clf
subplot(2,1,1)
plot(n,Y(1,1:L)-y0(1:L))
axis tight
xlabel('Sample number')
ylabel('Lowpass difference')
title('Difference in outputs produced by polyphase and direct forms')
subplot(2,1,2)
plot(n,Y(2,1:L)-y1(1:L))
axis tight
xlabel('Sample number')
ylabel('Highpass difference')
pause

% Plot the determinant of the polyphase matrix as a function of frequency.
R = 32;
W = 2/R*(-R/2:R/2-1);
H0even = fftshift(fft(H(1,1,:),R));
H0odd = fftshift(fft(H(1,2,:),R));
H1even = fftshift(fft(H(2,1,:),R));
H1odd = fftshift(fft(H(2,2,:),R));
delta = zeros(1,R);
delta(:) = H0even .* H1odd - H0odd .* H1even;
clf
plot(W,abs(delta),'x-')
axis([-1 1 0 1.5])
xlabel('Angular frequency (normalized by pi)')
ylabel('Magnitude of determinant')
title('Determinant of the polyphase matrix')
pause

% Verify that the filter is orthogonal i.e. Hp'(w*) Hp(w) = I
A11 = zeros(1,R);
A11(:) = abs(H0even).^2 + abs(H1even).^2;
A12 = zeros(1,R);
A12(:) = conj(H0even).*H0odd + conj(H1even).*H1odd;
A21 = zeros(1,R);
A21(:) = conj(H0odd).*H0even + conj(H1odd).*H1even;
A22 = zeros(1,R);
A22(:) = abs(H0odd).^2 + abs(H1odd).^2;
clf
subplot(4,1,1)
plot(W,A11,'x-');
axis([-1 1 0 1.5])
ylabel('A11')
title('Variation of A = Hp''(w*) Hp(w) with frequency')
subplot(4,1,2)
plot(W,A12,'x-')
axis([-1 1 0 1.5])
ylabel('A12')
subplot(4,1,3)
plot(W,A21,'x-')
axis([-1 1 0 1.5])
ylabel('A21')
subplot(4,1,4)
plot(W,A22,'x-')
axis([-1 1 0 1.5])
xlabel('Angular frequency (normalized by pi)')
ylabel('A22')



