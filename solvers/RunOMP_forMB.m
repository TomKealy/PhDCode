% OMP alogrithm for MMV/SMV: Y=AX which comes from the CTF (for real-valued MB signal)
% This function returns a solution Xr together with the support (Supp)
% The termination criteria are:
%    Stop if more than N bands were identified (N is the total number of bands in both parts of the spectrum range)
%    Stop if residual norm is greater than ResThreshold
%    Stop if the ratio between residual norm and the norm of the current solution is greater than ResvsSolThreshold
% Each iteration picks one support item. The symmetric location is added
% automatically. Once the support contains 2N locations, we check  adjacent
% locations. Finally we repeat "MoreIters" (with symmetric completion as
% before).
function [Supp] = RunOMP_forMB(Y,A,N, MoreIters)
[mA,nA] = size(A);
[mY,nY] = size(Y);

residual = Y;
OneSidedSupp = [];
Supp = [];
SymmetricSupp = [];
iter = 1;

NormACols = sqrt(diag(A'*A));

while (iter <= N+MoreIters)
    if  (iter <= N)
        % Matching step
        Z_1 = A'*residual;
        Z = sqrt(sum(abs(Z_1).^2,2))./NormACols;
        [maxVal, maxPos] = max(Z);
        BestLoc = maxPos(1);
        % Update support
        OneSidedSupp = [OneSidedSupp BestLoc];
        SymmetricLoc =  (nA+1-BestLoc);
        if (BestLoc ~= SymmetricLoc)
            SymmetricSupp = [SymmetricSupp SymmetricLoc];
        end
        Supp = [OneSidedSupp  SymmetricSupp];
    end
    if (iter == N)          %  Now we search for the adjacent positions
        CurrentLen = length(OneSidedSupp);
        for i=1:CurrentLen
            DominantBand = OneSidedSupp(i);
            CheckSupp = [DominantBand-1   DominantBand+1];
            CheckSupp = min(CheckSupp,nA);
            CheckSupp = max(CheckSupp,1);
            CheckSupp = unique(CheckSupp);
            if  (   (not(isempty(intersect(CheckSupp,OneSidedSupp))))   || ...
                    (not(isempty(intersect(CheckSupp,SymmetricSupp)))) )
                % meaning that we already picked the adjacent position
                continue;
            end
            % Matching step
            Z_2 = A(:,CheckSupp)'*residual;
            Z = sqrt(sum(abs(Z_2).^2,2))./NormACols(CheckSupp);
            [maxVal, maxPos] = max(Z);
            BestLoc = CheckSupp(maxPos(1));
            % Update support
            OneSidedSupp = [OneSidedSupp BestLoc];
            SymmetricLoc =  (nA+1-BestLoc);
            if (BestLoc ~= SymmetricLoc)
                SymmetricSupp = [SymmetricSupp SymmetricLoc];
            end
            Supp = [OneSidedSupp  SymmetricSupp];
        end
    end
    % Project residual
    As = A(:,Supp);
    solution = As*pinv(As)*Y;
    residual = Y-solution;

    iter = iter+1;
end

% Do not construct solution, since we need only the support
%Xr = zeros(nA,nY);
%Xr(Supp,:) = pinv(As)*Y;