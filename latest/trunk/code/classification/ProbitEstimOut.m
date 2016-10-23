% CLASS: ProbitEstimOut
% 
% HIERARCHY (Enumeration of the various super- and subclasses)
%   Superclasses: EstimOut
%   Subclasses: N/A
% 
% TYPE (Abstract or Concrete)
%   Concrete
%
% DESCRIPTION (High-level overview of the class)
%   The ProbitEstimOut class defines a scalar observation channel, p(y|z),
%   that constitutes a probit binary classification model, i.e., y is an
%   element of the set {0,1}, z is a real number, and
%             	 p(y = 1 | z) = Phi((z - Mean)/sqrt(Var)),
%   where Phi((x-b)/sqrt(c)) is the cumulative density function (CDF) of a
%   Gaussian random variable with mean b, variance c, and argument x.
%   Typically, mean = 0 and var = 1, thus p(y = 1 | z) = Phi(z).
%
% PROPERTIES (State variables)
%   Y           An M-by-T array of binary ({0,1}) class labels for the
%               training data, where M is the number of training data
%               points, and T is the number of classifiers being learned
%               (typically T = 1)
%   Mean        [Optional] An M-by-T array of probit function means (see
%               DESCRIPTION) [Default: 0]
%   Var         [Optional] An M-by-T array of probit function variances
%               (see DESCRIPTION).  Note that smaller variances equate to
%               more step-like sigmoid functions [Default: 1e-2]
%   maxSumVal 	Perform MMSE estimation (false) or MAP estimation (true)?
%               [Default: false]
%
% METHODS (Subroutines/functions)
%   ProbitEstimOut(Y)
%       - Default constructor.  Assigns Mean and Var to default values.
%   ProbitEstimOut(Y, Mean)
%       - Optional constructor.  Sets both Y and Mean.
%   ProbitEstimOut(Y, Mean, Var)
%       - Optional constructor.  Sets Y, Mean, and Var.
%   ProbitEstimOut(Y, Mean, Var, maxSumVal)
%       - Optional constructor.  Sets Y, Mean, Var, and maxSumVal.
%   estim(obj, zhat, zvar)
%       - Provides the posterior mean and variance of a variable z when
%         p(y|z) is the probit model and maxSumVal = false (see 
%         DESCRIPTION), and p(z) = Normal(zhat,zvar).  When maxSumVal =
%         true, estim returns MAP estimates of each element of z, as well
%         as the second derivative of log p(y|z).
%

%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 08/27/12
% Change summary: 
%       - Created (05/20/12; JAZ)
%       - Modified estim method to compute quantities using erfcx function
%         (v0.2) (07/10/12; JAZ)
%       - Added maxSumVal property to distinguish between MMSE and MAP
%         computations (08/27/12; JAZ)
% Version 0.2
%

classdef ProbitEstimOut < EstimOut
    
    properties
        Y;              % M-by-T vector of binary class labels
        Mean = 0;       % M-by-T vector of probit function means [dflt: 0]
        Var = 1e-2;    	% M-by-T vector of probit function variances [dflt: 1e-2]
        maxSumVal = false;   % Sum-product (false) or max-sum (true) GAMP?
    end
    
    methods
        % *****************************************************************
        %                      CONSTRUCTOR METHOD
        % *****************************************************************
        function obj = ProbitEstimOut(Y, Mean, Var, maxsumval)
            obj = obj@EstimOut;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.Y = Y;      % Set Y
                if nargin >= 2 && ~isempty(Mean)
                    % Mean property is an argument
                    obj.Mean = Mean;
                end
                if nargin >= 3 && ~isempty(Var)
                    % Var property is an argument
                    obj.Var = Var;
                end
                if nargin >= 4 && ~isempty(maxsumval)
                    % maxSumVal property is an argument
                    obj.maxSumVal = maxsumval;
                end
                if any(Var(:) < 0) && ~obj.maxSumVal
                    error('Var must be non-negative when running sum-product GAMP')
                elseif any(Var(:) <= 0) && obj.maxSumVal
                    error('Var must be strictly positive when running max-sum GAMP')
                end
            end
        end
        
        
        % *****************************************************************
        %                           SET METHODS
        % *****************************************************************
        function obj = set.Y(obj, Y)
            if ~all((Y(:) == -1) | (Y(:) == 0) | (Y(:) == 1))
                warning(['Elements of Y must be either in {0,1}' ...
                    ' or {-1,1}.  If you pass the true Z, then only use' ...
                    ' the genRand function or errors will occur.']);
            elseif any(Y(:) == -1)
                Y(Y == -1) = 0;
            end
            obj.Y = Y;
        end
        
        function obj = set.Mean(obj, Mean)
                obj.Mean = double(Mean);
        end
        
        function obj = set.Var(obj, Var)
            if any(Var(:) < 0)
                error('Var must be non-negative')
            else
                obj.Var = double(Var);
            end
        end
        
        function obj = set.maxSumVal(obj, maxsumval)
            if isscalar(maxsumval) && islogical(maxsumval)
                obj.maxSumVal = maxsumval;
            else
                error('ProbitEstimOut: maxSumVal must be a logical scalar')
            end
        end
        
        
        % *****************************************************************
        %                          ESTIM METHOD
        % *****************************************************************
        % This function will compute the posterior mean and variance of a
        % random vector Z whose prior distribution is N(Phat, Pvar), given
        % observations Y obtained through the separable channel model:
        % p(Y(m,t) = 1 | Z(m,t)) = Phi((Z(m,t) - Mean(m,t))/sqrt(Var(m,t)))
        % if obj.maxSumVal = false, otherwise it will return Zhat = argmax
        % log p(y|z) - 1/2/Pvar (z - Phat)^2 and the second derivative of
        % log p(y|z) evaluated at Zhat in Zvar, if obj.maxSumVal = true
        function [Zhat, Zvar] = estim(obj, Phat, Pvar)
            switch obj.maxSumVal
                case false
                    % Return sum-product expressions to GAMP in Zhat and
                    % Zvar
                    
                    % Start by computing the critical constant, C, on which
                    % the remainder of the computations depend.  Modulate 
                    % this constant by -1 for cases where Y(m) = 0.
                    PMonesY = sign(obj.Y - 0.1);	% +/- 1 for Y(m,t)'s
                    C = PMonesY .* ((Phat - obj.Mean) ./ ...
                        sqrt(Pvar + obj.Var));
                    
                    % Now compute the ratio normpdf(C)/normcdf(C)
                    ratio = (2/sqrt(2*pi)) * (erfcx(-C / sqrt(2)).^(-1));
                    
                    % Finally, compute E[Z(m,t) | Y(m,t)] = Zhat, and
                    % var{Z(m,t) | Y(m,t)} = Zvar
                    Zhat = Phat + PMonesY .* ...
                        ((Pvar ./ sqrt(Pvar + obj.Var)) .* ratio);
                    Zvar = Pvar - ((Pvar.^2 ./ (Pvar + obj.Var)) .* ...
                        ratio) .* (C + ratio);
                    
                case true
                    % Return max-sum expressions to GAMP in Zhat and Zvar
                    PMonesY = sign(obj.Y - 0.1);	% +/- 1 for Y(m,t)'s
                    
                    % Determine the expansion point about which to perform
                    % the Taylor series approximation
                    EP = (sign(PMonesY) == sign(Phat)) .* Phat;
                    
%                     % First compute a second-order Taylor series
%                     % approximation of log p(y|z) - 1/2/Pvar (z - Phat)^2,
%                     % about the point EP, and set as Zhat the maximizer 
%                     % of this approximation
%                     C = PMonesY .* (EP - obj.Mean) ./ sqrt(obj.Var);
%                     % Now compute the ratio normpdf(C)/normcdf(C)
%                     ratio = (2/sqrt(2*pi)) * (erfcx(-C / sqrt(2)).^(-1));
%                     % Compute 1st deriv of maximization functional
%                     Deriv1 = (PMonesY ./ sqrt(obj.Var)) .* ratio - ...
%                         (1./Pvar) .* (EP - Phat);
%                     % Compute 2nd deriv of maximization functional
%                     Deriv2 = -(1./obj.Var) .* ratio .* (C + ratio) - ...
%                         (1 ./ Pvar);
%                     % Set maximizer of Taylor approximation as Zhat
%                     Zhat = Phat - Deriv1 ./ Deriv2;

                    % Manually locate the value of z that sets the cost
                    % function derivative to zero using fsolve
                    opts = optimset('Jacobian', 'on', 'MaxIter', 25, ...
                        'Display', 'off');
                    F = @(z) zero_deriv(obj, z, Phat, Pvar);
                    Zhat = fsolve(F, EP, opts);
                    
                    % Now compute second derivative of log p(y|z) evaluated
                    % at Zhat (Note: Not an approximation)
                    % Start by computing frequently appearing constant
                    C = PMonesY .* (Zhat - obj.Mean) ./ sqrt(obj.Var);
                    % Now compute the ratio normpdf(C)/normcdf(C)
                    ratio = (2/sqrt(2*pi)) * (erfcx(-C / sqrt(2)).^(-1));
                    % Compute 2nd deriv of log p(y|z)
                    Deriv = -(1./obj.Var) .* ratio .* (C + ratio);
%                     Deriv = max(1e-6, Deriv);
                    
                    % Output in Zvar a function of the 2nd derivative that,
                    % once manipulated by gampEst, yields the desired
                    % max-sum expression for -g'_{out}
                    Zvar = Pvar ./ (1 - Pvar.*Deriv);
            end
        end
        
        
        % *****************************************************************
        %                         LOGLIKE METHOD
        % *****************************************************************
        % This function will compute *an approximation* to the expected
        % log-likelihood, E_z[log p(y|z)] when performing sum-product GAMP
        % (obj.maxSumVal = false).  The approximation is based on Jensen's 
        % inequality, i.e., computing log E_z[p(y|z)] instead.  If
        % performing max-sum GAMP (obj.maxSumVal = true), logLike returns
        % log p(y|z) evaluated at z = Zhat
        function ll = logLike(obj, Zhat, Zvar)
            PMonesY = sign(obj.Y - 0.1);	% +/- 1 for Y(m,t)'s
            switch obj.maxSumVal
                case false
                    % Start by computing the critical constant, C, on which
                    % the remainder of the computations depend.  Modulate 
                    % this constant by -1 for cases where Y(m,t) = 0.
                    C = PMonesY .* ((Zhat - obj.Mean) ./ sqrt(Zvar + obj.Var));
                    CDF = normcdf(C,0,1);
                    ll = log(CDF);
                    %Find bad values that cause log cdf to go to infinity
                    I = find(C < -30);
                    %This expression is equivalent to log(normpdf(C)) for
                    %all negative arguments greater than -38 and is
                    %numerically robust for all values smaller.  DO NOT USE
                    %FOR LARGE positive x.
                    ll(I) = -log(2)-0.5*C(I).^2+log(erfcx(-C(I)/sqrt(2)));
                case true
                    %Compute true constant
                    C = PMonesY .* (Zhat - obj.Mean)/sqrt(obj.Var);
                    ll = log(normcdf(C, 0, 1));
                    %Find bad values that cause log cdf to go to infinity
                    I = find(C < -30);
                    %This expression is equivalent to log(normpdf(C)) for
                    %all negative arguments greater than -38 and is
                    %numerically robust for all values smaller.  DO NOT USE
                    %FOR LARGE positive x.
                    ll(I) = -log(2)-0.5*C(I).^2+log(erfcx(-C(I)/sqrt(2)));
            end
        end
        
        % Compute output cost:
        % For sum-product compute
        %   (Axhat-phatfix)^2/(2*pvar) + log int_z p_{Y|Z}(y|z) N(z;phatfix, pvar) 
        %   with phatfix such that Axhat=estim(phatfix,pvar).
        % For max-sum GAMP, compute
        %   log p_{Y|Z}(y|z) @ z = Axhat
        function ll = logScale(obj,Axhat,pvar,phat)
                   
%           error('logScale method not implemented for this class. Set the GAMP option adaptStepBethe = false.');  

            %Find sign needed to compute the inner factor
            PMonesY = sign(obj.Y - 0.1);	% +/- 1 for Y(m,t)'s

            if~(obj.maxSumVal)
                % Find the fixed-point of phat
                opt.phat0 = Axhat;
                opt.alg = 1; % approximate newton's method
                opt.maxIter = 40; 
                opt.tol = 1e-4; 
                opt.stepsize = 0.2; 
                %Smallest regularization setting without getting the warnings for Var= 0
                opt.regularization = 1e-6; 
                opt.debug = false;
                phatfix = estimInvert(obj,Axhat,pvar,opt);     
                
                C = PMonesY .* ((phatfix - obj.Mean) ./ sqrt(pvar + obj.Var));
                % Compute log int_z p_{Y|Z}(y|z) N(z;phatfix, pvar)
                % Has a closed form solution: see C. E. Rasmussen, 
                % "Gaussian Processes for Machine Learning." sec 3.9.
                ls = log(normcdf(C, 0 , 1));
                %Find bad values that cause log cdf to go to infinity
                I = find(C < -30);
                %This expression is equivalent to log(normpdf(C)) for
                %all negative arguments greater than -38 and is
                %numerically robust for all values smaller.  DO NOT USE
                %FOR LARGE positive x.
                ls(I) = -log(2)-0.5*C(I).^2+log(erfcx(-C(I)/sqrt(2)));
                
                % Combine to form output cost
                ll = ls + 0.5*(Axhat - phatfix).^2./pvar;
            else
                %Compute true constant
                C = PMonesY .* (Axhat - obj.Mean)/sqrt(obj.Var);
                ll = log(normcdf(C, 0, 1));
                %Find bad values that cause log cdf to go to infinity
                I = find(C < -30);
                %This expression is equivalent to log(normpdf(C)) for
                %all negative arguments greater than -38 and is
                %numerically robust for all values smaller.  DO NOT USE
                %FOR LARGE positive x.
                ll(I) = -log(2)-0.5*C(I).^2+log(erfcx(-C(I)/sqrt(2)));
            end
            
        end
        
        
        % *****************************************************************
        %                       NUMCOLUMNS METHOD
        % *****************************************************************
        function S = numColumns(obj)
            % Return number of columns of Y
            S = size(obj.Y, 2);
        end
        
        %Generate random samples given the distribution
        function y = genRand(obj, z)
            if obj.Var > 0
                actProb = normcdf(z,obj.Mean, sqrt(obj.Var));
                y = (rand(size(z)) < actProb);
            elseif obj.Var == 0
                y = sign(z-obj.Mean);
                y(y == -1) = 0;
            else
                error('Var must be non-negative')
            end
        end
    end
    
    methods (Access = private)
        % *****************************************************************
        %                         FSOLVE METHOD
        % *****************************************************************
        % This method is used by MATLAB's fsolve function in the max-sum
        % case to locate the value of z that sets the prox-operator
        % derivative to zero
        function [Fval, Jacobian] = zero_deriv(obj, z, phat, pvar)
            % Compute value of derivative at z
            PMonesY = sign(obj.Y - 0.1);	% +/- 1 for Y(m,t)'s
            C = PMonesY .* (z - obj.Mean) ./ sqrt(obj.Var);
            % Now compute the ratio normpdf(C)/normcdf(C)
            ratio = (2/sqrt(2*pi)) * (erfcx(-C / sqrt(2)).^(-1));
            % Value of derivative
            Fval = PMonesY.*ratio./sqrt(obj.Var) - (z - phat)./pvar;
            
            % Optionally compute Jacobian of F at z
            if nargout >= 2
                M = numel(phat);
                Jvec = -(1./obj.Var) .* ratio .* (C + ratio) - (1 ./ pvar);
%                 Jacobian = diag(Jvec);      % Okay for small problems
                Jacobian = sparse(1:M, 1:M, Jvec);
            end
        end
    end
end