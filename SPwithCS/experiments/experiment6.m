%% experiment6.m - Produces data for Figure 9a
%
% Most recent change - 9/15/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

function experiment6(fileno)

filename = ['experiment6_data' num2str(fileno)];

addpath('../');

% Set signal/algorithm parameters
N = 1024*4;
B = 1/256;
K = 5;

dpssSizeMin = N*B; % Initial value of dpssSize - should be less than N*B, but >=1
dpssSizeMax = 38;
dpssWidth = N*B/2;  % DPSS width
        
kkmax = 50; % Number of iterations per parameter setting to test
mmmax = 13;
M_max = N*B*K*mmmax;
jjmax = mmmax*3;

%senseOpts.sampler = 'matrix';
%randmat = randn(M_max,N);
senseOpts.sampler = 'rdmdemod';
senseOpts.pnseq = sign(randn(N,1)); % Choose random pn sequence

sigOpts.N = N; sigOpts.B = B; sigOpts.K = K;
sigOpts.type = 'randtones';  % Random tones model
sigOpts.TPB = 50;  % 50 tones per band

cosampOpts.N = N;
cosampOpts.maxiter = 20;
cosampOpts.tol = 1e-3;
cosampOpts.modulator = generateModulator(N,B); % Precompute modulator matrix
cosampOpts.verbose = 0;
cosampOpts.detectAlg = 'blockOMP';
cosampOpts.debias = 1;

signalNorm = zeros(kkmax,1);
errorNorm1 = zeros(mmmax,kkmax);
errorNorm2 = inf*ones(mmmax,kkmax,jjmax);

timevec = zeros(mmmax,kkmax);

for kk=1:kkmax,
    
    [x supp] = generateSignal(sigOpts); % Generate signal
    signalNorm(kk) = norm(x); % Store signal norm
    cosampOpts.normBound = signalNorm(kk)+1e-13;
        
    for mm=1:mmmax,
        
        tic
        
        M = N*B*K*(mm+1)/2;
        senseOpts.subsamp = N/M;
        
        %  DPSS KNOBS
        dpssSize = round(5.5*M/(N*B*K)+5);
        if(dpssSize > dpssSizeMax)
            dpssSize = dpssSizeMax;
        end
        if(dpssSize < dpssSizeMin)
            dpssSize = dpssSizeMin;
        end
        cosampOpts.PsiRecover = dpss(N,dpssWidth,dpssSize); % Precompute DPSS basis
        
%        senseOpts.matrix = 1/sqrt(M)*randmat(1:M,:);
        y = sense(x,0,senseOpts); % Compute measurements
        
        [x_hat, ~, ~] = dpssCosamp(y,senseOpts,K,B,cosampOpts);
        errorNorm1(mm,kk) = norm(x - x_hat);
        
        for jj=mm:3*mm,
            [x_fft,~] = fftCosamp(y,senseOpts,K*jj,cosampOpts);
            errorNorm2(mm,kk,jj) = norm(x - x_fft);
        end
        
        timevec(mm,kk) = toc;
        
        [minerr,minval] = min(errorNorm2(mm,kk,:));
        disp(['mm = ' num2str(mm) ' SNR: ' num2str(20*log10(signalNorm(kk)/minerr)) ' K selected: '  num2str(K*minval)  ' --  Time: ' num2str(timevec(mm,kk)/60) ' min'])
       
    end
    disp(['Elapsed time for iteration ' num2str(kk) ' : ' num2str(sum(timevec(:,kk))/60) ' min'])
    save(filename,'signalNorm','errorNorm1','errorNorm2','timevec');
end
