function params = updateParameters(rhat, vr, initParams, estfun, varargin)
% UPDATEPARAMETERS computes parameters for AWGN estimation problem.
%
% params = updateParameters(rhat, vr, initParams, estfun, varargin)
%
% Input:
% - rhat: noisy measurements rhat = x + w, where w_i~N(0, vr) (n x 1)
% - vr: variance of AWGN (n x 1)
% - initParams: initial value of parameters to estimate
% - estfun: cell array of handlers to scalar function used for estimation
% - varargin: option for the algorithm
%       Options:
%       - nIterations: number of iterations (def: 50)
%       - plotReconstruction: detailed plots (def: false)
%       - trueParams: oracle parameters
%       - tol: stopping criterion (def: 1e-4)
%       - tolCount: maximum number of iterations to allow no change (def: 3)
%       - verbose: print to command line (def: false)
%
% Output:
% - params: updated parameters
%
% U. S. Kamilov, BIG, EPFL, 2012.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

absDiff = @(x, xprev) (norm(x(:)-xprev(:))^2)/(norm(xprev(:))^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of parameters
nParameters = length(initParams);

% Problem dimension
n = length(rhat);

% Number of iterations
nIterations = 50;

% Plot intermediate results
plotReconstruction = false;

% Display intermediate messages
verbose = false;

% Answer difference tolerance
tol = 1e-4;

% Number of allowed iteration within tol
tolCount = 3;

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
        case 'niterations'
            nIterations = value;
        case 'plotreconstruction'
            plotReconstruction = value;
        case 'trueparams'
            trueParams = value;
        case 'tol'
            tol = value;
        case 'tolcount'
            tolCount = value;
        case 'verbose'
            verbose = value;
        otherwise
            error('updateParameters: input is not recognized!');
    end
end

% Control parameters
assert(nIterations >= 0, 'UPDATEPARAMETERS: nIterations >= 0.');
assert(tol >= 0, 'UPDATEPARAMETERS: tol >= 0.');
assert(tolCount >= 0, 'UPDATEPARAMETERS: tolCount >= 0.');
assert(length(estfun) == nParameters,...
    'UPDATEPARAMETERS: length(estfun) == nParameters.');

if(exist('trueParams', 'var'))
    assert(length(trueParams) == nParameters,...
        'UPDATEPARAMETERS: length(trueParams) == nParameters.');
end

if(plotReconstruction && ~exist('trueParams', 'var'))
    plotReconstruction = false;
    fprintf('UPDATEPARAMETERS: trueParams not provided.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(plotReconstruction)
    h = figure('Color', 'w');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterative estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = initParams;
paramsprev = 1e6*ones(1, nParameters);

for iIteration = 1:nIterations
    
    % Go through each parameter
    for iParameter = 1:nParameters
        
        % Update function for this parameter
        fi = estfun{iParameter};
        
        % Update current parameter
        params(iParameter) = sum(fi(rhat, vr, params))/n;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Feedback
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % CMD
        if(verbose)
            fprintf('[Update Parameters: %d/%d]', iIteration, nIterations);
            if(exist('trueParams', 'var'))
                err = absDiff(params(iParameter)-trueParams(iParameter));
                fprintf('[p%d error: %.4f]', iParameter, err);
            end
            fprintf('\n');
        end
        
        % Plot
        if(plotReconstruction)
            err = absDiff(params(iParameter)-trueParams(iParameter));
            
            figure(h);
            subplot(nParameters, 1, iParameter);
            plot(iIteration, params(iParameter), 'b.');
            hold on;
            semilogy([1, nIterations],...
                [trueParams(iParameter), trueParams(iParameter)],...
                'r-', 'LineWidth', 1.2)
            title(sprintf('Parameter %d: %.4f. Error: %.4f',...
                iParameter, params(iParameter), err));
            xlim([1 nIterations]);
            drawnow;
        end
    end
    
    % If without a change decrese tolerance
    if(absDiff(params, paramsprev) < tol)
        tolCount = tolCount - 1;
    end
    
    % Stopping criterion
    if(tolCount <= 0)
        break;
    end
    
    % Update previous estimate
    paramsprev = params;
end