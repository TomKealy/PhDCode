function [eigspect,centerfreq_hat,x_hat,Rz,Phi_Omega]=mc_recovery(y,L,W,c,t,nyqover)
% This script recovers the band support and the signal amplitudes of the original multiband
% signal x from the multi-coset samples y. The recovery algorithm is that of Feng and Bresler
% (see the technical report "Multi-coset Sampling and Recover of Sparse Multiband 
% Signals", accompanying the CTSS Sampling Toolbox).

% NOTE: In this script, the size of the support (K_hat) is set manually. Its proper setting
% depends on the sparsity of the vector s(omega)---again, see the technical report 
% "Multi-coset Sampling and Recover of Sparse Multiband Signals".

% Usage: [eigspect,centerfreq_hat,x_hat,Rz,Phi_Omega]=mc_recovery(y,L,W,c,t,nyqover)
% y: output samples of multi-coset sampler
% L: parameter that determines the period of the nonuniform sampling (every L/W seconds q samples are collected)
% W: Nyquist rate (Hz)
% c: multi-coset sampling pattern
% K: number of center frequencies to be extracted from MUSIC spectrum
% t: vector of discrete time instances over which band sparse signal is defined
% eigspect: MUSIC eigenspectrum 
% centerfreq_hat: estimate of center frequencies of occupied bands (Hz)
% x_hat: reconstructed (estimated) signal

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

%% Support recovery (active frequency band recovery)  
q=size(y,2);  
yinterp=zeros(q,size(y,1)*L*nyqover); % initialize interpolated output
z=zeros(q,size(y,1)*L*nyqover); % initialize interpolated and delayed output

for i=1:q
  yinterp(i,:)=interpft(y(:,i),size(y,1)*L*nyqover); % interpolate by L
  z(i,c(i)*nyqover:size(y,1)*L*nyqover) = (L/W)*yinterp(i,1:size(y,1)*L*nyqover-(c(i)*nyqover)+1); % delay by sampling pattern c(i)
end

Rz=2*pi*(L/W)*(z*z'); % form covariance matrix
[eigvectors,eigvalues]=eig(Rz); % eigen-decomposition of R
[unused,i]=sort(diag(eigvalues),'ascend'); % sort eigenvalues in ascending order
%*******************
K_hat=12; % for now, choose K_hat manually 
%*******************
Un=eigvectors(:,i(1:q-K_hat)); % declared noise-only eigenvectors of R

m=ceil(-(L+1)/2)+1:floor((L+1)/2);
Phi=exp(-1i*(2*pi/L)*(c-1)*m); % form Phi matrix

eigspect=zeros(1,L); % initialize eigen-spectrum (MUSIC spectrum)
for i=1:L
 eigspect(i)=(Phi(:,i)'*(Un*Un')*Phi(:,i))/(Phi(:,i)'*Phi(:,i)); % form eigen-spectrum (project columns of Phi onto noise-only subspace)
end
[unused,bandindices]=sort(abs(eigspect),'ascend');
Phi_Omega=Phi(:,bandindices(1:K_hat)); % extract the columns of Phi that lie in the null space of the noise-only subspace 

% Convert bandindices into corresponding center frequencies (Hz)
centfreq=(-W/2:W/L:W/2);
centerfreq_hat(:,1)=centfreq(L-bandindices(1:K_hat)+1);

%% Reconstruct simulated analog signals
s=(W/L)*pinv(Phi_Omega)*z;
x_hat=sum(s.*exp(1i*2*pi*centerfreq_hat*t(1:size(s,2))),1);