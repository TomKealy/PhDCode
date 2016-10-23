function y = emVarx(rhat, vr, params)
% EMVARX scalar function used for EM estimation of varx
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
assert(rho >= 0 && rho <= 1, 'emVarx: rho > 0 && rho < 1');
assert(varx >= 0, 'emVarx: varx >= 0');

svar = vr+varx;

gr = 1./(1 + ((1/rho) - 1)*(sqrt(1 + varx./vr)).*exp(-0.5*(rhat.^2).*(1./vr - 1./(vr+varx))));

y = ((((rhat .* varx ./ svar).^2) + (varx .* vr ./ svar)) .* gr)/rho;