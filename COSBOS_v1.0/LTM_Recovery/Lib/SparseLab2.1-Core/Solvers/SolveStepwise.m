function [betahat, activationHist, tHist] = SolveStepwiseNoStoppingCriterion(X,y,threnter,intercept)
% SolveStepwise: Forward Stepwise Regression
% Usage
%       [beta activationHist] = SolveStepwise(X,y,threnter)
% Input
%
% Outputs
%    betahat          solution of Stepwise
%    activationHist   Array of indices showing elements entering the solution set
%    thrHist          Array of t-statistics entering the solution set
% Description
%   SolveStepwise is a forward stepwise regression algorithm
%   which starts without model terms and adds in terms a greedy fashion.
%

if nargin < 3
    threnter = 4;
end

if nargin < 4
    intercept = 1;
end

[n,p] = size(X);
if intercept
    ActiveSet = ones(n,1); % intercept vector
else
    ActiveSet = [];
end
ActiveIndices = [1:p];
activationHist = [];
tHist = [];
betahat = zeros(p,1);
done = 0;

while ~done

    thr = -Inf;

    % for each available index
    for i = ActiveIndices

        % least squares
        Xnew = [ActiveSet X(:,i)];
        b = Xnew \ y;

        % standard error
        ehat = y - Xnew*b;
        s2 = sum(ehat.^2);
        if (n ~= length(b))
            s2 = s2./(n-length(b));
        end
        se = sqrt(diag(inv(Xnew'*Xnew))*s2);

        % t-statistic
        if (se ~= 0)
            t = b./se;
            if abs(t(end)) > thr
                ind = i;
                thr = t(end);
            end
        end
    end

    % criterion to enter the model
    if (abs(thr) >= threnter)
        activationHist = [activationHist ind];
        tHist = [tHist thr];
        ActiveSet = [ActiveSet X(:,ind)];
        ActiveIndices = setdiff(ActiveIndices,ind);
    else
        done = 1;
    end

end


%
% Copyright (c) 2006. Victoria Stodden
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

