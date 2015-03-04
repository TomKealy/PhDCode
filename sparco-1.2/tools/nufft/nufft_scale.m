 function sn = nufft_scale(Nd, Kd, alpha, beta)
%function sn = nufft_scale(Nd, Kd, alpha, beta)
% Compute scaling factors for NUFFT
% in:
%	Nd,Kd
%	alpha	{d}
%	beta	{d}
% out:
%	sn	[[Nd]]		scaling factors
%
% Copyright 2004-7-8, Jeff Fessler, The University of Michigan

if nargin < 4, help(mfilename), nufft_scale_test, error(mfilename), end

dd = length(Nd);
if dd == 1 & ~iscell(alpha) % 1D case
	sn = nufft_scale1(Nd(1), Kd(1), alpha, beta);
	return;
end

%
% scaling factors: "outer product" of 1D vectors
%
sn = 1;
for id=1:dd
	tmp = nufft_scale1(Nd(id), Kd(id), alpha{id}, beta{id});
	sn = sn(:) * tmp';
end
if length(Nd) > 1
	sn = reshape(sn, Nd);	% [(Nd)]
else
	sn = sn(:);	% [(Nd)]
end


% Compute scaling factors for 1D NUFFT (from Fourier series coefficients)
% in:
%	N,K
%	alpha
%	beta
% out:
%	sn	[N]		scaling factors
%
% Copyright 2001-10-4, Jeff Fessler, The University of Michigan

function sn = nufft_scale1(N, K, alpha, beta)

if ~isreal(alpha(1)), error 'need real alpha_0', end
L = length(alpha) - 1;

%
% compute scaling factors from Fourier coefficients
%
if L > 0
	sn = zeros(N,1);
	n = [0:(N-1)]';
	i_gam_n_n0 = i * (2*pi/K) * (n - (N-1)/2) * beta;

	for l1=-L:L
		alf = alpha(abs(l1)+1);
		if l1 < 0, alf = conj(alf); end
		sn = sn + alf * exp(i_gam_n_n0 * l1);
	end

else
	sn = alpha * ones(N,1);
end

if 0
	printf('range real(sn) = %g,%g', min(real(sn)), max(real(sn)))
	printf('range imag(sn) = %g,%g', min(imag(sn)), max(imag(sn)))
end


%
% self test
%
function nufft_scale_test

N = 100;
K = 2*N;
alpha = [1.0 -0.0 -0.2];
sn = nufft_scale(N, K, alpha, 1);
clf, plot(1:N, real(sn), 'y-', 1:N, imag(sn), 'g-')
legend('sn real', 'sn imag')
