%NPR_SYNTHESIS Near perfect reconstruction synthesis
%  Y = NPR_SYNTHESIS(C,X) synthesizes a wideband signal Y from a
%  number of subbands stored in X. Each subband is a row in X.
%  C is a two dimensional array, containing the filter coefficients.
%  The number of rows in X must be twice the number of rows in C. 
%
% (c) 2007 Wessel Lubberhuizen

function [y] = npr_synthesis(coeff,x)

% number of channels
N=size(coeff,1);

% number of slices
M=size(x,2);

% split into even and odd channels
x = reshape(x,2,N,M);

y1 = squeeze(x(1,:,:));
y2 = squeeze(x(2,:,:));

% apply dft
y1 = fft(y1,[],1)*N;
y2 = fft(y2,[],1)*N;

% apply channel filters
for i=1:N
    y1(i,:) = filter(coeff(i,:),1,y1(i,:));
    y2(i,:) = filter(coeff(i,:),1,y2(i,:));
end

% apply frequency shift
for i=1:N;
    y2(i,:) = y2(i,:) * exp(-sqrt(-1)*pi*(i-1)/N);
    y2(i,2:2:M) = -y2(i,2:2:M);
end

% combine filter results
y = reshape(y1-y2,1,M*N);
