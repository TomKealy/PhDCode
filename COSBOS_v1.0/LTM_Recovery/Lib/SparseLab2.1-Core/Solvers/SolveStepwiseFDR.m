function [betahat, activationHist, tHist] = SolveStepwiseFDR(X,y,intercept,q)
% SolveStepwiseFDR: Forward Stepwise Regression using FDR stopping rule
% Usage
%       [beta activationHist] = SolveStepwiseFDR(X,y,intercept)
% Input
%       q             FDR parameter 0<q<1
% Outputs
%    betahat          solution of Stepwise
%    activationHist   Array of indices showing elements entering the solution set
%    thrHist          Array of t-statistics entering the solution set
% Description
%   SolveStepwise is a forward stepwise regression algorithm
%   which starts without model terms and adds in terms a greedy fashion.
%

if nargin < 3
    intercept = 1;
end

if nargin < 4
    q=.25;
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

    % FDR statistic
        k = size(activationHist,2)+1; % number of terms in actual model
        FDRvalue = q*k/p; % p is number of potential terms in the model
        
        df = n - k;
        
        tlargest = 0;
        ilargest = 0;
        
    % for each available index
    for i = ActiveIndices

        % least squares
        Xnew = [ActiveSet X(:,i)];
        Xnew;
        stats = regstats(y,Xnew,eye(k),'tstat');
        %bhat = Xnew \ y %cannot use QR because n < p
        bhat = stats.tstat.beta;
        
        t = stats.tstat.t(end);
        if (isnan(t) || t==0), 
            t=0;
            pval=1; 
        end
        if abs(t) > tlargest,
            if (t==Inf), 
                tlargest = 1000000;
                ilargest = i;
                pval=0; 
            else
                tlargest = abs(t);
                ilargest = i;
                pval = stats.tstat.pval(end);
            end
        end;
    end

    % criterion to enter the model
    % p-value

    if (pval < FDRvalue)
        activationHist = [activationHist ilargest];
        tHist = [tHist tlargest];
        ActiveSet = [ActiveSet X(:,ilargest)];
        ActiveIndices = setdiff(ActiveIndices,ilargest);
    else
        done = 1;
    end
end
ActiveSet;

if intercept
    Xj = [ones(n,1) X(:,activationHist)];
    betahattilda = Xj \ y;
    betahat(activationHist) = betahattilda(2:end);
else
    betahattilda = X(:,activationHist) \ y;
    betahat(activationHist) = betahattilda;
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

