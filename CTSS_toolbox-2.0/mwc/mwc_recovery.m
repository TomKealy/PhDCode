function [centerfreq_hat,x_hat,R,Psi,A]=mwc_recovery(y,L,W,Phi,M,t)
% This script recovers the band support and the signal amplitudes of the original multiband
% signal x from samples y acquired by the Modulated Wideband Converter (MWC). The recovery 
% algorithm is that of Mishali and Eldar (see the technical report "Sampling Sparse Multiband 
% Signals with a Modulated Wideband Converter", accompanying the CTSS Sampling Toolbox).

% Usage: [centerfreq_hat,x_hat,R,Psi,A]=mwc_recovery(y,L,W,Phi,M,t)
% y: output samples of MWC
% L: specifies rate and period of random +/- 1 sequence wrt Nyquist rate (1/Tp=W/L)
% W: Nyquist rate (Hz)
% Phi: random +/- 1 matrix
% M: specifies sampling rate per channel (1/Ts=W/M)
% t: vector of discrete time instances over which band sparse signal is defined
% centerfreq_hat: estimate of center frequencies of active bands (Hz)
% x_hat: reconstructed signal (signal estimate)

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

%% Support recovery (active frequency band recovery) 
R=2*pi*(M/W)*(y*y'); % estimate covariance matrix
[eigvectors,eigvalues]=eig(R); % eigen-decomposition of R

l(:,1)=0:L-1;    m=ceil((-L/(2*M))*(M+1))+1:floor((L/(2*M))*(M+1));
Psi=exp(-1i*(2*pi/L)*l*m); % form Psi matrix
alpha=ones(size(m))*1/L; % complex factors do not appear in digital simulation (see technical 
                         % report "Sampling Sparse Multiband Signals with a Modulated Wideband Converter")
A=Phi*Psi*diag(alpha); % matrix characterising linear system of equations

V=eigvectors*sqrt(eigvalues);
[S,Supp,ResNorm,NormResvsSol]=RunOMP_Unnormalized(V,A,100,0,0.01,'false'); % find support by solving MMV via OMP (M. Mishali's algorithm)
K_hat=length(Supp);  
A_Omega=A(:,Supp(1:K_hat)); % reduce dimension of A 

% Associate support indices to corresponding center frequencies (Hz)
centfreq= -W/2:W/L:W/2;
centerfreq_hat(:,1)=centfreq(Supp(1:K_hat)+1);

%% Reconstruct Nyquist samples                            
s_hat = pinv(A_Omega)*y; % y and s_hat are matrices of time domain signals
s_interp = zeros(size(s_hat,1),length(t));
for i=1:size(s_interp,1) 
 s_interp(i,:) = interpft(s_hat(i,:),length(t)); % interpolate to expand to original (Nyquist) length
end
 
x_hat=sum(s_interp.*exp(-1i*2*pi*centerfreq_hat*t),1); % modulate spectrum slices to appropriate locations