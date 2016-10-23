% This demp compares MSE predicted by SE analysis with the MSE
% computed by oracle and adaptive GAMP. We also plot LASSO.
%
% Corresponding paper:
% U. S. Kamilov, S. Rangan, A. K. Fletcher, and M. Unser, 
% "Approximate Message Passing with Consistent Parameter Estimation and 
% Applications to Sparse Learning," 
% IEEE Trans. Inf. Theory., vol. 60, no. 5, pp. 2969-2985, May 2014.
%
% U. S. Kamilov, BIG, EPFL, 2012.
% Web-site: http://bigwww.epfl.ch/kamilov

clear; close all; home;

addpath(genpath(sprintf('.%s', filesep)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Length of the signal
n = 400;

% Sparsity ratio
rho = 0.2;

% Signal variance
varx = 5;

% Measurement ratios (beta = m/n)
nBetas = 4;
minBeta = 0.5;
maxBeta = 2;
betas = linspace(minBeta, maxBeta, nBetas);

% Noise variances
varn = 0.1;

% Number of problem realizations
nInstances = 100;

% Number of iteration of the algorithm
nIterations = 50;
nEmIterations = 100;

% Initialization
xhat0 = 0*ones(n, 1);
vx0 = rho*varx*ones(n, 1);
rho0 = 0.5;
varx0 = 50;

% Adaptive initialization
init0 = {xhat0, vx0, [rho0, varx0]};

% Non-adaptive initialization
init = {xhat0, vx0, [rho, varx]};

% Estimation function
estfun{1} = @ (rhat, vr, params) emRho(rhat, vr, params);
estfun{2} = @ (rhat, vr, params) emVarx(rhat, vr, params);
denfun = @(rhat, vr, params) awgnBgShrink(rhat, vr, params);

% MSE functions
computeMse = @(noise) 10*log10((norm(noise(:))^2)/length(noise));
computeSnr = @(sig, noise) 10*log10((norm(sig(:))^2)/(norm(noise(:))^2));

% LASSO params
lmin = 1e-3;
lmax = 100;

mse1 = zeros(nInstances, nBetas); % oracle GAMP
mse2 = zeros(nInstances, nBetas); % adaptive GAMP
mse3 = zeros(nInstances, nBetas); % LASSO

time1 = zeros(nInstances, nBetas); % oracle GAMP
time2 = zeros(nInstances, nBetas); % adaptive GAMP
time3 = zeros(nInstances, nBetas); % LASSO


for iInstance = 1:nInstances
    for iBeta = 1:nBetas
        
        % Number of measurements
        m = round(n*betas(iBeta));
        
        % Beta
        beta = m/n;
        betas(iBeta) = beta;

        % Generate the signal (Gaussian with proba. rho)
        x = binornd(1, rho, n, 1) .* randn(n, 1) .* sqrt(varx);
        
        % Generate the measurement matrix
        A = (1/sqrt(m)) .* randn(m, n);
        
        % Obtain measurements
        z = A*x;
        
        % Noisy observations
        y = z + sqrt(varn) * randn(m, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run GAMP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Oracle
        trueParams = {x, [rho, varx]};
        
        fprintf('[%d/%d][%d/%d] GAMP: ', iBeta, nBetas, iInstance, nInstances);
        
        t = cputime;
        
        % Run non-adaptive GAMP
        [xhat1, outs1] = adaptiveGamp(A, y, init, denfun, estfun,...
            'isAdaptive', false,...
            'nIterations', nIterations,...
            'nParamIterations', nEmIterations,...
            'NoiseVariance', varn,...
            'oracle', trueParams,...
            'Tol', 1e-4);
        
        % Store results
        time1(iInstance, iBeta) = cputime-t;
        mse1(iInstance, iBeta) = outs1.mse(end);        
        
        fprintf('[MSE = %.4f][T = %.4f]\n', outs1.mse(end), time1(iInstance, iBeta));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run Adaptive GAMP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('[%d/%d][%d/%d] Adaptive GAMP: ', iBeta, nBetas, iInstance, nInstances);
        
        t = cputime;
        
        % Run adaptive GAMP
        [xhat2, outs2] = adaptiveGamp(A, y, init0, denfun, estfun,...
            'isAdaptive', true,...
            'nIterations', nIterations,...
            'nParamIterations', nEmIterations,...
            'NoiseVariance', varn,...
            'oracle', trueParams,...
            'Tol', 1e-4);
        
        % Store results
        time2(iInstance, iBeta) = cputime-t;
        mse2(iInstance, iBeta) = outs2.mse(end);        
        
        fprintf('[MSE = %.4f][T = %.4f]\n', outs2.mse(end), time2(iInstance, iBeta));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run LASSO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf('[%d/%d][%d/%d] LASSO: ', iBeta, nBetas, iInstance, nInstances);
        
        % Run LASSO
        [~, lam] = oracle(@(dat,l) l1_ls(A, dat,l,1e-3, 1), x, y, lmin, lmax);
        
        t = cputime;
        
        xhat3 = l1_ls(A, y, lam,1e-3, 1);
        
        % Store results
        time3(iInstance, iBeta) = cputime-t;
        mse3(iInstance, iBeta) = computeMse(xhat3-x);        
        
        fprintf('[MSE = %.4f][T = %.4f]\n', mse3(iInstance, iBeta), time3(iInstance, iBeta));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w', 'Name', 'Measurements Rate vs. MSE');
plot(betas, mean(mse3, 1), 'dm-',...
    betas, mean(mse1, 1), 'sr--',...
    betas, mean(mse2, 1), 'ob-', 'LineWidth', 1.2);
set(gca,'FontSize',14);
xlabel('Measurement ratio (m/n)', 'fontsize', 14);
ylabel('MSE (dB)', 'fontsize', 14);
legend('Lasso', 'Oracle', 'Adaptive', 'Location', 'NorthEast');
xlim([betas(1), betas(end)]);