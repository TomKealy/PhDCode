function [x supp] = generateSignal(opts)
%% Generates a multiband signal according to the options specified
% 
% Usage: [x supp] = generateSignal(opts)
% Inputs: opts - a struct setting the signal parameters
%           opts.type - 'dpss', 'randtones', 'randproc'
%           opts.N - signal length
%           opts.K - number of active bands
%           opts.B - bandwidth of each band 
%                    (B must be <= .5 and 1/B must be an integer)
%           opts.supp - optional parameter to specify the support pattern 
%           ***other options available depending on signal type***
% Outputs: x - generated multiband signal
%          supp - vector indicating occupied frequency bands
%
% Most recent change - 8/27/2011
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

%% Initialize parameters
N = opts.N; B = opts.B; K = opts.K;

% Set support based on opts or generate it randomly
if(isfield(opts,'supp'))
    supp = opts.supp;
else
    idx = randperm(1/B);
    supp = idx(1:K);
end

x = zeros(N,1);  % Initialize x

%% Generate something that will live within our DPSS model by construction
if strcmp(opts.type,'dpss')
    % Generate modulator matrix
    modulator = generateModulator(N,B);
    % Construct base DPSS basis
    if(isfield(opts,'dpssWidth'))
        dpssWidth = opts.dpssWidth;
    else
        dpssWidth = N*B/2;
    end
    if(isfield(opts,'dpssSize')) 
        dpssSize = opts.dpssSize;
    else
        dpssSize = 2*dpssWidth;
    end
    Psi = dpss(N,dpssWidth,dpssSize);
    % Generate test signal supported on supp 
    basis = repmat(Psi,1,K).*kron(modulator(:,supp),ones(1,dpssSize));
    alpha = randn(size(basis,2),1) + 1j*randn(size(basis,2),1);
    x = basis*alpha;
end

%% Generate a random sum of tones uniformly distributed in the bands
if strcmp(opts.type,'randtones')
    t_vec = [0:1:N-1];  % Construct time index
    % Select number of tones per band
    if(isfield(opts,'TPB'))
        TPB = opts.TPB;
    else
        TPB = 50;
    end
    for jj=1:length(supp),
        random_freqs = ((supp(jj)-1)*B -.5) + B*rand(1,TPB);
        tones = exp(2j*pi*t_vec'*random_freqs);
        x = x + tones*(randn(TPB,1) + 1j*randn(TPB,1));
    end
end

%% Generate a signal as a random process which is white on the bands
if strcmp(opts.type,'randproc')
    % Oversampling factor to use in generating longer signal to truncate
    if(isfield(opts,'oversample'))
        OS = opts.oversample;
    else
        OS = 8;
    end
    x = zeros(OS*N,1);
    for jj=1:length(supp),
        noise = randn(OS*N,1) + 1i*randn(OS*N,1);
        NOISE = fft(noise);
        FILT = zeros(OS*N,1);
        f1 = (supp(jj)-1)*B-.5;  % Define active band
        f2 = f1+B;
        if(f1 < 0)  % Translate [-.5,.5] to [0,1]
            k1 = OS*N*(f1+1) +1;
            if(f2 < 0)
                k2 = OS*N*(f2+1) +1;
                FILT(k1:k2) = 1;
            else
                k2 = OS*N*f2+1;
                FILT(1:k2) = 1;
                FILT(k1:end) = 1;
            end
        else
            k1 = OS*N*f1 + 1;
            k2 = OS*N*f2 + 1;
            FILT(k1:k2) = 1;
        end
        x = x + ifft(NOISE.*FILT);
    end
    x = x(OS*N/2:OS*N/2+N-1);
end