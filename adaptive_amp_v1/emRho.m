function y = emRho(rhat, vr, params)
% EMRHO scalar function used for EM estimation of rho.
%
% Input:
% - rhat: noisy measurements rhat = x + w, w~N(0, vr) (n x 1)
% - vr: noise variance (n x 1)
% - params: initial parameters [rho, varx]
%
% Output:
% - y: array (n x 1) containing pointwise estimates
%
% U. S. Kamilov, BIG, EPFL, 2012.

% Extract parameters
rho = params(1);
varx = params(2);

% Check validity
assert(rho >= 0 && rho <= 1, 'emRho: rho > 0 && rho < 1');
assert(varx >= 0, 'emRho: varx >= 0');

y = 1./(1 + ((1/rho) - 1)*(sqrt(1 + varx./vr)).*exp(-0.5*(rhat.^2).*(1./vr - 1./(vr+varx))));