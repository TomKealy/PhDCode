function [xhat, vx] = awgnBgShrink(rhat, vr, params)
%BGSHRINK returns MMSE estimate from AWGN measurements of a GaussBernoulli
%random variable.
%
% [xhat, vx] = bgShrink(rhat, vr, rho, varx)
%
% Rhat = X + V, where V ~ N(0, vr) and X ~ rho*N(0, varx) + (1-rho)*dirac(x)
%
% Input:
% - rhat: Noisy measurement (n x 1)
% - vr: AWGN variance (n x 1)
% - params: [rho, variance], where rho is sparsity (e.g. 0.1) and varx is
%       variance of the prior
%
% Output:
% - xhat = E[x | Rhat = rhat]
% - vx = var [x | Rhat = rhat]
%
% U. S. Kamilov, BIG, EPFL, 2012.

rho = params(1);
varx = params(2);

assert(rho <= 1 && rho >= 0, 'rho <= 1 && rho >= 0');
assert(all(vr) >= 0, 'vr >= 0');
assert(varx >= 0, 'varx >= 0');

% Posterior support proba
if(rho >= 1)
    gr = 1;
elseif(rho <= 0)
    gr = 0;    
else
    gr = 1./(1 + ((1/rho) - 1)*(sqrt(1 + varx./vr)).*exp(-0.5*(rhat.^2).*(1./vr - 1./(vr+varx))));
end

xhat = gr .* rhat .* varx ./ (vr + varx);
vx = gr .* ((varx .* vr ./ (varx + vr)) + (rhat .* varx ./ (vr + varx)).^2)-...
    (gr .* rhat .* varx ./ (vr + varx)).^2;