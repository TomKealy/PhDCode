function [Xr, Supp, ResNorm, NormResvsSol] = RunOMP_Unnormalized(Y,A,NumIters, ResThreshold, ResvsSolThreshold, SymmetricSupport)
% OMP alogrithm for MMV/SMV: Y=AX
% This function returns a solution Xr together with the support (Supp)
% The termination criteria are:
%    Stop if more than NumIters iterations
%    Stop if residual norm is greater than ResThreshold
%    Stop if the ratio between residual norm and the norm of the current solution is greater than ResvsSolThreshold
% If SymmetricSupport is true, then the algorithm select pairs of indices every
%      iteration, such that the indices are symettric with respect to solution
%       dimensions (special feature for conjugate symmetric solutions).
[mA,nA] = size(A);
[mY,nY] = size(Y);

residual = Y;
Supp = [];
iter = 1;
resnorm = norm(Y);
resvssolnorm = inf;

NormACols = sqrt(diag(A'*A));

while (  (iter <= NumIters)  &&  (resnorm > ResThreshold) && (resvssolnorm >  ResvsSolThreshold))
    % Matching step
    Z_1 = A'*residual;
    Z = sqrt(sum(abs(Z_1).^2,2))./NormACols;
    [maxVal, maxPos] = max(Z);
    BestLoc = maxPos(1);
    % Update support
    Supp = [Supp BestLoc];
    if (SymmetricSupport == true)
        SymmetricLoc =  (nA+1-BestLoc);
        if (BestLoc ~= SymmetricLoc)
            Supp = [Supp SymmetricLoc];
        end
    end
    % Project residual
    As = A(:,Supp);
    solution = As*pinv(As)*Y;
    residual = Y-solution;
    % update norms
    resnorm = norm(residual);
    resvssolnorm  = resnorm/norm(solution);
    ResNorm(iter) = resnorm;
    NormResvsSol(iter) = resvssolnorm ;
    % increment iter number
    iter = iter+1;
end

% Construct solution
Xr = zeros(nA,nY);
Xr(Supp,:) = pinv(As)*Y;