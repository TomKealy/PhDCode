%% Demonstration of Signal Space CoSaMP with non-normalized orthogonal dictionary
%
% Most recent change - 7/31/2012
%
% Copyright 2012, Mark Davenport, Deanna Needell, Michael Wakin
%
% This file is part of Signal Space CoSaMP Toolbox version 1.0.
%
%    Signal Space CoSaMP Toolbox is free software: you can redistribute it 
%    and/or modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    Signal Space CoSaMP Toolbox is distributed in the hope that it will be 
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Signal Space CoSaMP Toolbox.  If not, 
%    see <http://www.gnu.org/licenses/>.
%
% Running this software requires the following packages to be pre-installed:
%
%    CVX, available at http://cvxr.com/cvx/download/
%
% For exporting a pdf plot [optional], the following m-files are needed:
%
%    plotpdftex.m, available at http://www.mathworks.de/matlabcentral/fileexchange/20847
%    laprint.m, available at http://www.mathworks.de/matlabcentral/fileexchange/4638



%% Initialize
clear all; close all;

%% Set signal/algorithm parameters
n = 256;
mVec = [20:8:100];
k = 8;

% Dictionary: Orthobasis with non-normalized columns
d = n;
diagVec = [100*ones(n/2,1); ones(n/2,1)];
D = diag(diagVec);

numTrials = 1000;

clear errors alphas
errors{1} = zeros(numTrials,length(mVec));
errors{2} = zeros(numTrials,length(mVec));
errors{3} = zeros(numTrials,length(mVec));
alphas{1} = zeros(numTrials,length(mVec),d);
alphas{2} = zeros(numTrials,length(mVec),d);
alphas{3} = zeros(numTrials,length(mVec),d);
alphas{4} = zeros(numTrials,length(mVec),d); % true alpha

for trialNum = 1:numTrials
    for mIndex = 1:length(mVec)
        m = mVec(mIndex);
        disp(['Beginning trial ' num2str(trialNum) ', m = ' num2str(m)]);
        
        % Choose random measurement matrix
        A = randn(m,n)/sqrt(m);
        
        % Generate test signal
        alpha = zeros(d,1);
        idx = randperm(d);
        alpha(idx(1:k)) = randn(k,1);
        x = D*alpha;
        
        % Generate measurements
        y = A*x; % Compute measurements
        noise = randn(m,1);
        noise = noise/norm(noise);
        yn = y; % + 10^(-4)*noise;
        
        % Signal Space Cosamp Opts
        cosampOpts.maxiter = 20;
        cosampOpts.tol = 1e-6;
        cosampOpts.verbose = 0;
        cosampOpts.normBound = 1e6*norm(alpha);
        cosampOpts.alg = 'nnd';
        
        % Single matrix opts
        singleMatrixOpts.normBound = 1e6*norm(alpha);
        
        %% Reconstruct
        
        % Signal Space CoSaMP with optimal projections
        [x1,alpha1] = doubleMatrixCoSaMP(yn,A,D,k,cosampOpts);
        
        % Traditional CoSaMP
        alpha2 = singleMatrixCoSaMP(yn,A*D,k,singleMatrixOpts);
        x2 = D*alpha2;
        
        % Traditional OMP
        alpha3 = singleMatrixOMP(yn,A*D,k,singleMatrixOpts);
        x3 = D*alpha3;
        
        %% Display results
        disp([' Double-matrix Cosamp Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x1))) 'dB']);
        disp([' Cosamp Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x2))) 'dB']);
        disp([' OMP Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x3))) 'dB']);
        
        errors{1}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x1));
        errors{2}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x2));
        errors{3}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x3));
        alphas{1}(trialNum,mIndex,:) = alpha1;
        alphas{2}(trialNum,mIndex,:) = alpha2;
        alphas{3}(trialNum,mIndex,:) = alpha3;
        alphas{4}(trialNum,mIndex,:) = alpha;
        
    end
end

% save demo_nnd_results

if 0
    load demo_nnd_results
        
    numTrials = trialNum; % in case the saved file was an intermediate one from the simulation
    medianThreshold = round(0.5*numTrials);
    ninetyFivePercentThreshold = round(0.05*numTrials);
    perfectRecoveryThreshold = 100; % in dB
    clear medianErrors ninetyFivePercentErrors pctRecoveryBySNR pctRecoveryBySupport
    for testType = 1:3
        medianErrors{testType} = zeros(length(mVec),1);
        ninetyFivePercentErrors{testType} = zeros(length(mVec),1);
        pctRecoveryBySNR{testType} = zeros(length(mVec),1);
        for mIndex = 1:length(mVec)
            eTemp{testType} = sort(squeeze(errors{testType}(1:numTrials,mIndex)),'ascend');
            medianErrors{testType}(mIndex) = eTemp{testType}(medianThreshold);
            ninetyFivePercentErrors{testType}(mIndex) = eTemp{testType}(ninetyFivePercentThreshold);
            pctRecoveryBySNR{testType}(mIndex) = sum(eTemp{testType} > perfectRecoveryThreshold)/numTrials*100;
            correctSupps = 0;
            for ii = 1:numTrials                
                estSupp = find(alphas{testType}(ii,mIndex,:));
                trueSupp = find(alphas{4}(ii,mIndex,:));
                if isequal(estSupp,trueSupp)
                   correctSupps = correctSupps + 1;
                end                
            end
            pctRecoveryBySupport{testType}(mIndex) = correctSupps/numTrials*100;
        end        
    end

    figure(1); clf; hold on;
    kIndex = 1;
    plot(mVec, pctRecoveryBySNR{1}, 'kx-', 'LineWidth', 2);
    plot(mVec, pctRecoveryBySNR{2}, 'r+-', 'LineWidth', 2);
    xlabel('$m$','FontSize',20);
    ylabel('Percent of perfect recovery','FontSize',20);
    h = legend('Signal Space CoSaMP','CoSaMP','Location','SouthEast');
    set(h,'FontSize',24);
    set(gca,'FontSize',20);
    plotpdftex(1,'demo_nnd.pdf', [1.5 1])
end