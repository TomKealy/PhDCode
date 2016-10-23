function [xhat, outs] = adaptiveGamp(A, y, init, denfun, estfun, varargin)
% ADAPTIVEGAMP GAMP combined with parameter estimation.
% Solves AWGN problem: y = A*x + w, where w_i~N(0, varw).
%
% [xhat, outs] = adaptiveGamp(A, y, init, denfun, estfun, varargin)
%
% Input:
% - A: measurement matrix (m x n)
% - y: noisy measurements: y = A*x + w (m x 1)
% - init: initialization for the signal and parameters: {xhat0, vx0,
%       theta0}, where theta0 contains parameters of the prior
% - denfun: prior dependent denoising function (handle)
% - estfun: array of function handles to update parameters theta of the
%       prior.
% - varargin: options for the algorithm
%       Options:
%       - isAdaptive: tell whether to adapt (def: true)
%       - nIterations: number of iterations for GAMP (def: 50)
%       - nParamIterations: number of EM iterations (def: 500)
%       - NoiseVariance: AWGN variance (def: 0)
%       - oracle: original signal and parameters {x, theta}.
%       - PlotReconstruction: GAMP reconstruction plots (def: false)
%       - PlotParamReconstruction: update plots (def: false)
%       - Tol: stopping criterion for the algorithms (def: 1e-4)
%       - TolCount: max. num. times Tol must be satisfied to stop (def: 3)
%       - TolParam: stopping criterion for parameter update (def: 1e-4)
%       - TolParamCount: max. num. times for TolParam to stop (def: 3)
%       - verbose: print command line messages (def: false)
%
% Output:
% - xhat: final reconstruction (n x 1)
% - outs: extra data
%       - outs.mse: per iteration MSE (requires oracle)
%       - outs.time: per iteration CPU-time
%       - outs.params: estimated parameters
%
% U. S. Kamilov, BIG, EPFL, 2012.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computeMse = @(noise) 10*log10((norm(noise(:))^2)/length(noise));
computeSnr = @(sig, noise) 10*log10((norm(sig(:))^2)/(norm(noise(:))^2));
absDiff = @(x) sum(abs(x(:)))/length(x(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of parameters
nParameters = length(init{3});

% Number of GAMP iterations
nIterations = 50;

% Number of EM iterations
nParamIterations = 500;

% AWGN variance
noiseVariance = 0;

% Tolerances
Tol = 1e-4;
TolParam = 1e-4;

% Tolerance counter
TolCount = 3;
TolParamCount = 3;

% Plot the progress of GAMP reconstruction
plotReconstruction = false;

% Plot parameter estimation reconstructions progress
plotParamReconstruction = false;

% Display intermediate messages
verbose = false;

% Is oracle set?
isOracleSet = false;

% Is adaptive?
isAdaptive = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of options
nargs = length(varargin);

% Go through options
for i = 1:2:nargs
    % Extract name/value pair
    name = lower(varargin{i});
    value = varargin{i+1};
    
    switch(name)
        case 'tolcount'
            TolCount = value;
        case 'tolparamcount'
            TolParamCount = value;
        case 'niterations'
            nIterations = value;
        case 'nparamiterations'
            nParamIterations = value;
        case 'noisevariance'
            noiseVariance = value;
        case 'oracle'
            oracle = value;
            isOracleSet = true;
        case 'isadaptive'
            isAdaptive = value;
        case 'plotreconstruction'
            plotReconstruction = value;
        case 'plotparamreconstruction'
            plotParamReconstruction = value;
        case 'tol'
            Tol = value;
        case 'tolparam'
            TolParam = value;
        case 'verbose'
            verbose = value;
        otherwise
            error('adaptiveGamp: input is not recognized!');
    end
end

% Control parameters
assert(noiseVariance >= 0, 'AdGAMP: noiseVariance >= 0.');
assert(Tol >= 0, 'AdGAMP: Tol >= 0.');
assert(TolCount >= 0, 'AdGAMP: TolCount >= 0.');
assert(TolParam >= 0, 'AdGAMP: TolParam >= 0.');
assert(TolParamCount >= 0, 'AdGAMP: TolParamCount >= 0.');
assert(nIterations >= 0, 'AdGAMP: nIterations >= 0.');
assert(nParamIterations >= 0, 'AdGAMP: nParamIterations >= 0.');
assert(length(estfun) == nParameters,...
    'AdGAMP: length(estfun) == nParameters.');

if(~isOracleSet)
    plotReconstruction = false;
    plotParamReconstruction = false;
    fprintf('AdGAMP: oracle not provided. Cannot plot.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem dimensions
[m, n] = size(A);

% Hadamard product of the matrix
AA = A.*A;

% Previous estimate
xhatprev = 1e6*ones(n, 1);

% Initialize GAMP
shat = zeros(m, 1);
xhat = init{1};
vx = init{2};

% Track CPU time per iteration
outs.time = zeros(nIterations, 1);

% Initialize EM
params = init{3};

if(isOracleSet)
    % Extract oracle
    x = oracle{1};
    trueParams = oracle{2};
    
    outs.mse = zeros(nIterations, 1);
    outs.snr = zeros(nIterations, 1);
end

if(plotReconstruction)
    h = figure('Color', 'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iIteration = 1:nIterations
    
    % Message to command line
    if(verbose)
        fprintf('--------------------\n');
        fprintf('[GAMP: %d/%d]\n', iIteration, nIterations);
        fprintf('--------------------\n');
    end
    
    % Track time
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Factor update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Linear
    vp = AA*vx;
    phat = A*xhat - vp.*shat;
    
    % Non-Linear
    vs = 1 ./ (vp+noiseVariance);
    shat = vs .* (y - phat);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Variable update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Linear
    vr = 1./(AA' * vs);
    rhat = xhat + vr .* (A' * shat);
    
    if(isAdaptive)
        % Estimate parameters
        if(isOracleSet)
            params = updateParameters(rhat, vr, params, estfun,...
                'nIterations', nParamIterations,...
                'trueParams', trueParams,...
                'tol', TolParam,...
                'tolCount', TolParamCount,...
                'verbose', verbose,...
                'plotReconstruction', plotParamReconstruction);
        else
            params = updateParameters(rhat, vr, params, estfun,...
                'nIterations', nParamIterations,...
                'tol', TolParam,...
                'tolCount', TolParamCount,...
                'verbose', verbose,...
                'plotReconstruction', plotParamReconstruction);
        end
    end
    
    % Non-Linear
    [xhat, vx] = denfun(rhat, vr, params);
    
    % Store time and parameters
    outs.time(iIteration) = toc;
    outs.params = params;
    
    % Compute error metrics
    if(isOracleSet)
        outs.mse(iIteration) = computeMse(xhat - x);
        outs.snr(iIteration) = computeSnr(x, xhat-x);
    end
    
    % CMD
    if(verbose)
        if(isOracleSet)
            fprintf('[MSE: %.4f][SNR: %.4f]\n',...
                outs.mse(iIteration), outs.snr(iIteration));
        end
        fprintf('\n');
    end
    
    % Plot
    if(plotReconstruction)
        figure(h);
        
        subplot(3, 1, 1);
        plot(iIteration, outs.mse(iIteration), 'b.');
        hold on;
        title(sprintf('MSE(%d/%d): %.f dB', iIteration, nIterations, outs.mse(iIteration)));
        xlim([1 nIterations]);
        drawnow;
        
        subplot(3, 1, 2);
        plot(iIteration, outs.snr(iIteration), 'b.');
        hold on;
        title(sprintf('SNR(%d/%d): %.f dB', iIteration, nIterations, outs.snr(iIteration)));
        xlim([1 nIterations]);
        drawnow;
        
        subplot(3, 1, 3);
        plot(1:n, x, 'ro', 1:n, xhat, 'b.');
        title('Estimate');
        xlim([1 n]);
        drawnow;
    end
    
    % If without a change decrese tolerance
    if(absDiff(xhat - xhatprev) < Tol)
        TolCount = TolCount - 1;
    end
    
    % Stopping criterion
    if(TolCount <= 0)
        outs.snr = outs.snr(1:iIteration);
        outs.mse = outs.mse(1:iIteration);
        outs.time = outs.time(1:iIteration);
        break;
    end
    
    % Update previous estimate
    xhatprev = xhat;
end