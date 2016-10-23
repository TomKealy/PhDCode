classdef SoftThreshDMMEstimIn < EstimIn
    % SoftThreshDMMEstimIn:  Performs Donoho/Maleki/Montanari-style soft thresholding
    % when using GAMP to solve "min_x 1/2/var*norm(y-A*x,2)^2 + lambda*norm(x,1)". 
    % Warning: the val output is currently set to zero
    
    properties
        %DMM threshold gain "alpha", which is NOT the "lambda" in the cost function. 
	%Here, the soft threshold is set as thresh = alpha * sqrt(mean(mur));
        alpha;
    end
    
    methods
        % Constructor
        function obj = SoftThreshDMMEstimIn(alpha)
            obj = obj@EstimIn;
            if nargin ~= 0 % Allow nargin == 0 syntax
                obj.alpha = alpha;
            end
        end
        
        % Prior mean and variance
        function [mean0, var0, valInit] = estimInit(~)
            mean0 = 0; %For now, set these to arbitrary constants
            var0  = 1e-2;
            valInit = -Inf;
        end
        
        % Carry out soft thresholding
        function [xhat,xvar,val] = estim(obj,rhat,rvar)
            
            
            %Compute the threshold
            thresh = obj.alpha * sqrt(mean(rvar));
            
            %Estimate the signal
            xhat = max(0,abs(rhat)-thresh) .* sign(rhat);
            
            %Estimate the variance
            xvar = rvar .* (mean(double(abs(xhat) > 0))*ones(size(xhat)));
            
            %For now, let's set the val output to zero. Not clear
            %how to properly set this without knowledge of "lambda"
            val = zeros(size(xhat));
            
            
        end
        
    end
    
end

