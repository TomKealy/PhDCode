%% Demonstration of Signal Space CoSaMP with 4x overcomplete DFT dictionary and with sparse coefficient vectors having well-separated (not clustered) nonzeros
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
%    SPGL1, available at http://www.cs.ubc.ca/~mpf/spgl1/download.html
%
% For exporting a pdf plot [optional], the following m-files are needed:
%
%    plotpdftex.m, available at http://www.mathworks.de/matlabcentral/fileexchange/20847
%    laprint.m, available at http://www.mathworks.de/matlabcentral/fileexchange/4638


%% Initialize
clear all; close all;

%% Set signal/algorithm parameters
n = 256;
k = 8;
mVec = [20:8:100];

% Dictionary: Overcomplete DFT
d = 4*n;
D = ifft(eye(d))*d/sqrt(n);
D = D(1:n,:);

numTrials = 500;
clear errors alphas
errors{1} = zeros(numTrials,length(mVec));
errors{2} = zeros(numTrials,length(mVec));
errors{3} = zeros(numTrials,length(mVec));
errors{4} = zeros(numTrials,length(mVec));
errors{5} = zeros(numTrials,length(mVec));
errors{6} = zeros(numTrials,length(mVec));
alphas{1} = zeros(numTrials,length(mVec),d);
alphas{2} = zeros(numTrials,length(mVec),d);
alphas{3} = zeros(numTrials,length(mVec),d);
alphas{4} = zeros(numTrials,length(mVec),d);
alphas{5} = zeros(numTrials,length(mVec),d);
alphas{6} = zeros(numTrials,length(mVec),d);
alphas{7} = zeros(numTrials,length(mVec),d); % true alpha

for trialNum = 1:numTrials
    for mIndex = 1:length(mVec)
        m = mVec(mIndex);
        
        disp(['Beginning trial ' num2str(trialNum) ', m = ' num2str(m)]);
        
        % Choose random measurement matrix
        A = randn(m,n)/sqrt(m);
        
        % Generate test signal with well-separated nonzeros in coefficient
        % vector
        minSep = 0;
        while minSep < 9
            supp = randperm(d);
            supp = sort(supp(1:k),'ascend');
            minSep = min(diff(supp));
        end
        alpha = zeros(d,1);
        alpha(supp) = randn(k,1) + 1i*randn(k,1);
        x = D*alpha;
        
        % Generate measurements
        y = A*x; % Compute measurements
        noise = randn(m,1);
        noise = noise/norm(noise);
        yn = y; % + 10^(-4)*noise;
        
        % Signal Space Cosamp Opts
        doubleMatrixOpts.maxiter = 20;
        doubleMatrixOpts.tol = 1e-6;
        doubleMatrixOpts.verbose = 0;
        doubleMatrixOpts.normBound = 10*norm(alpha);
        
        % Single matrix opts
        singleMatrixOpts.maxiter = 20; % used only by singleMatrixCoSaMP
        singleMatrixOpts.normBound = 10*norm(alpha);
        
        %% Reconstruct
        
        % Signal Space CoSaMP with CoSaMP for computing projections
        doubleMatrixOpts.alg = 'cosamp';
        [x1,alpha1] = doubleMatrixCoSaMP(yn,A,D,k,doubleMatrixOpts);
        
        % Signal Space CoSaMP with OMP for computing projections
        doubleMatrixOpts.alg = 'omp';
        [x2,alpha2] = doubleMatrixCoSaMP(yn,A,D,k,doubleMatrixOpts);
        
        % Traditional CoSaMP
        alpha3 = singleMatrixCoSaMP(yn,A*D,k,singleMatrixOpts);
        x3 = D*alpha3;
        
        % Traditional OMP
        alpha4 = singleMatrixOMP(yn,A*D,k,singleMatrixOpts);
        x4 = D*alpha4;
        
        % Signal Space CoSaMP with l1-minimization for computing projections
        doubleMatrixOpts.alg = 'l1';
        [x5,alpha5] = doubleMatrixCoSaMP(yn,A,D,k,doubleMatrixOpts);
        
        % Traditional l1-minimization
        spgOpts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
        alpha6 = spg_bp(A*D, y, spgOpts);             % Solves alpha6 = argmin ||alpha||_1 s.t. A*D*alpha = y
        [alphaSort,alphaInds] = sort(abs(alpha6),'descend');
        supp6 = alphaInds(1:k);
        alpha6 = zeros(d,1);
        alpha6(supp6) = lsOnSupport(x,D(:,supp6),singleMatrixOpts.normBound); % Debiasing step after L1
        x6 = D*alpha6;
        
        disp([' Signal Space CoSaMP (CoSaMP) Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x1))) 'dB']);
        disp([' Signal Space CoSaMP (OMP) Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x2))) 'dB']);
        disp([' CoSaMP Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x3))) 'dB']);
        disp([' OMP Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x4))) 'dB']);
        disp([' Signal Space CoSaMP (L1) Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x5))) 'dB']);
        disp([' L1 Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x6))) 'dB']);
        
        errors{1}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x1));
        errors{2}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x2));
        errors{3}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x3));
        errors{4}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x4));
        errors{5}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x5));
        errors{6}(trialNum,mIndex) = 20*log10(norm(x)/norm(x-x6));
        alphas{1}(trialNum,mIndex,:) = alpha1;
        alphas{2}(trialNum,mIndex,:) = alpha2;
        alphas{3}(trialNum,mIndex,:) = alpha3;
        alphas{4}(trialNum,mIndex,:) = alpha4;
        alphas{5}(trialNum,mIndex,:) = alpha5;
        alphas{6}(trialNum,mIndex,:) = alpha6;
        alphas{7}(trialNum,mIndex,:) = alpha;
    end
    if mod(trialNum,25) == 0
        save demo_4xdft_separated_results
    end
end

% save demo_4xdft_separated_results

if 0
    load demo_4xdft_separated_results
    
    numTrials = trialNum; % in case the saved file was an intermediate one from the simulation        
    medianThreshold = round(0.5*numTrials);
    ninetyFivePercentThreshold = round(0.05*numTrials);
    perfectRecoveryThreshold = 100; % in dB
    clear medianErrors ninetyFivePercentErrors pctRecoveryBySNR pctRecoveryBySupport
    for testType = 1:6
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
                trueSupp = find(alphas{7}(ii,mIndex,:));
                if isequal(estSupp,trueSupp)
                   correctSupps = correctSupps + 1;
                end                
            end
            pctRecoveryBySupport{testType}(mIndex) = correctSupps/numTrials*100;
        end        
    end
    
    colorVec = {'k','b','r','g','c','m'};
    figure(1); clf;
    subplot(2,1,1); hold on;
    for testType = 1:6
        plot(mVec,medianErrors{testType}(:),char(colorVec(testType)));
    end
    xlabel('m');
    ylabel('SNR');
    title(['Median performance']);
    legend('SSCoSaMP (CoSaMP)','SSCoSaMP (OMP)','CoSaMP','OMP','SSCoSaMP (L1)','L1','Location','Best');
    subplot(2,1,2); hold on;
    for testType = 1:6
        plot(mVec,ninetyFivePercentErrors{testType}(:),char(colorVec(testType)));
    end
    xlabel('m');
    ylabel('SNR');
    title(['95 percent performance']);
    legend('SSCoSaMP (CoSaMP)','SSCoSaMP (OMP)','CoSaMP','OMP','SSCoSaMP (L1)','L1','Location','Best');
    
    figure(2); clf; hold on;
    for testType = 1:6
        plot(mVec,pctRecoveryBySNR{testType}(:),char(colorVec(testType)));
    end
    xlabel('m');
    ylabel('Pct recovery');
    title(['Percent of recovery above ' num2str(perfectRecoveryThreshold) 'dB SNR']);
    legend('SSCoSaMP (CoSaMP)','SSCoSaMP (OMP)','CoSaMP','OMP','SSCoSaMP (L1)','L1','Location','Best');
    
    figure(3); clf; hold on;
    for testType = 1:6
        plot(mVec,pctRecoveryBySupport{testType}(:),char(colorVec(testType)));
    end
    xlabel('m');
    ylabel('Pct recovery');
    title(['Percent of recovery with perfect support']);    
    legend('SSCoSaMP (CoSaMP)','SSCoSaMP (OMP)','CoSaMP','OMP','SSCoSaMP (L1)','L1','Location','Best');
    
    figure(4); clf; hold on;
    colorVec = {'kx-','bo-','r+-','g*-','cs-','md-'};
    for testType = 1:6
        plot(mVec,pctRecoveryBySNR{testType}(:),char(colorVec(testType)),'LineWidth',2);
    end
    xlabel('$m$','FontSize',20);
    ylabel('Percent of perfect recovery','FontSize',20);
    h = legend('SSCoSaMP (CoSaMP)','SSCoSaMP (OMP)','CoSaMP','OMP','SSCoSaMP ($\ell_1$)','$\ell_1$','Location','SouthEast');
    set(h,'FontSize',16);
    set(gca,'FontSize',20);
    plotpdftex(4,'demo_separated.pdf', [1.5 1]); % have to enlarge legend box by hand before printing
    
end
