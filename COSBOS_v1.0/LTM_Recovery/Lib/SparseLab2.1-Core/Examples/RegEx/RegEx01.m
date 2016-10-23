function RegEx01(X,y,b)
% Plots results for Forward Stepwise Algorithm against the true
% underlying model.

%Find sparse solutions, with true model unknown

[bhatstep, activationHist, tHist]=SolveStepwiseFDR(X,y,0);

% Plot Results

figure
plot(1:200,bhatstep,'-');
hold on
plot(1:200,b,'.r')
title('Forward Stepwise Regression Solution')
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
