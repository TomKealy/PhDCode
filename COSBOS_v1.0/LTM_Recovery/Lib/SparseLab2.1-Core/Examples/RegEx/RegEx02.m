function RegEx02(X,y,b)
% Plots results for Forward Stepwise Algorithm, with False Discovery Rate Threshold,
% against the true underlying model.

%Find sparse solutions, with true model unknown

[bhatFDR, activationHist, tHist]=SolveStepwiseFDR(X,y,0);

% Plot Results

figure;
plot(1:200,bhatFDR,'-');
hold on
plot(1:200,b,'.r')
title('Forward Stepwise (with False Discovery Rate Threshold) Regression Solution')
xlabel('Variable Number')
ylabel('Coefficient Value')
legend('Estimated Coefficient Values','True Coefficient Values')

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
