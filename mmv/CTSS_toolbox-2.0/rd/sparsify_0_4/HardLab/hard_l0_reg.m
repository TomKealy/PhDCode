function [s, err_mse, iter_time]=hard_l0_reg(x,A,m,lambda,varargin)
% hard_l0_reg: Hard thresholding algorithm that keeps values with 
% magnitude greater than lambda in each iteration. 
% Guaranteed to find local optimum of 
%               min_s || x  - Ps ||2  + lambda^2 || s ||_0
% Requires || P ||_2 < 1 or specification of step size to ensure convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage
%   [s, err_mse, iter_time]=hard_l0_reg(x,P,m,lambda,'option_name','option_value')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%   Mandatory:
%               x   	 Observation vector to be decomposed
%               P   	 Either:
%                       	1) An nxm matrix (n must be dimension of x)
%                       	2) A function handle (type "help function_format" 
%                          	   for more information)
%                          	   Also requires specification of P_trans option.
%                       	3) An object handle (type "help object_format" for 
%                          	   more information)
%               m   	 length of s 
%               lambda   threshold to be used
%
%   Possible additional options:
%   (specify as many as you want using 'option_name','option_value' pairs)
%   See below for explanation of options:
%__________________________________________________________________________
%   option_name    |     available option_values                | default
%--------------------------------------------------------------------------
%   stopTol        | number (see below)                         | 1e-6
%   P_trans        | function_handle (see below)                | 
%   maxIter        | positive integer (see below)               | n^2
%   verbose        | true, false                                | false
%   start_val      | vector of length m                         | zeros
%   step_size      | number [0 1]                               | zeros
%
%  stopping criteria used : (OldRMS-NewRMS)/RMS(x) < stopTol
%
%   stopTol: Value for stopping criterion.
%
%   P_trans: If P is a function handle, then P_trans has to be specified and 
%            must be a function handle. 
%
%   maxIter: Maximum number of allowed iterations.
%
%   verbose: Logical value to allow algorithm progress to be displayed.
%
%   start_val: Allows algorithms to start from partial solution.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%    s              Solution vector 
%    err_mse        Vector containing mse of approximation error for each 
%                   iteration
%    iter_time      Vector containing computation times for each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%   Implements the l0-regularised algorithm described in [1].
%   This algorithm takes a gradient step and then thresholds to only retain
%   elements with a magnitude above lambda.
%   
% References
%   [1] T. Blumensath and M.E. Davies, "Iterative Thresholding for Sparse 
%       Approximations", submitted, 2007
% See Also
%   hard_l0_Mterm
%
% Copyright (c) 2007 Thomas Blumensath
%
% The University of Edinburgh
% Email: thomas.blumensath@ed.ac.uk
% Comments and bug reports welcome
%
% This file is part of sparsity Version 0.1
% Created: April 2007
%
% Part of this toolbox was developed with the support of EPSRC Grant
% D000246/1
%
% Please read COPYRIGHT.m for terms and conditions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Default values and initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[n1 n2]=size(x);
if n2 == 1
    n=n1;
elseif n1 == 1
    x=x';
    n=n2;
else
   error('x must be a vector.');
end
    
sigsize     = x'*x/n;
oldERR      = sigsize;
err_mse     = [];
iter_time   = [];
STOPTOL     = 1e-5;
MAXITER     = n^2;
verbose     = false;
initial_given=0;
s_initial   = zeros(m,1);
mu=1;

if verbose
   display('Initialising...') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargout 
    case 3
        comp_err=true;
        comp_time=true;
    case 2 
        comp_err=true;
        comp_time=false;
    case 1
        comp_err=false;
        comp_time=false;
    case 0
        error('Please assign output variable.')        
    otherwise
        error('Too many output arguments specified')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Look through options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Put option into nice format
Options={};
OS=nargin-4;
c=1;
for i=1:OS
    if isa(varargin{i},'cell')
        CellSize=length(varargin{i});
        ThisCell=varargin{i};
        for j=1:CellSize
            Options{c}=ThisCell{j};
            c=c+1;
        end
    else
        Options{c}=varargin{i};
        c=c+1;
    end
end
OS=length(Options);
if rem(OS,2)
   error('Something is wrong with argument name and argument value pairs.') 
end
for i=1:2:OS
   switch Options{i}
        case {'stopTol'}
            if isa(Options{i+1},'numeric') ; STOPTOL     = Options{i+1};   
            else error('stopTol must be number. Exiting.'); end
        case {'P_trans'} 
            if isa(Options{i+1},'function_handle'); Pt = Options{i+1};   
            else error('P_trans must be function _handle. Exiting.'); end
        case {'maxIter'}
            if isa(Options{i+1},'numeric'); MAXITER     = Options{i+1};             
            else error('maxIter must be a number. Exiting.'); end
        case {'verbose'}
            if isa(Options{i+1},'logical'); verbose     = Options{i+1};   
            else error('verbose must be a logical. Exiting.'); end 
        case {'start_val'}
            if isa(Options{i+1},'numeric') && length(Options{i+1}) == m ;
                s_initial     = Options{i+1};   
                initial_given=1;
            else error('start_val must be a vector of length m. Exiting.'); end
        case {'step_size'}
            if isa(Options{i+1},'numeric') && (Options{i+1}) > 0 && (Options{i+1}) <=1 ;
                mu     = Options{i+1};   
            else error('Stepsize must be between 0 and 1. Exiting.'); end
        otherwise
            error('Unrecognised option. Exiting.') 
   end
end






if nargout >=2
    err_mse = zeros(MAXITER,1);
end
if nargout ==3
    iter_time = zeros(MAXITER,1);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Make P and Pt functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if          isa(A,'float')      P =@(z) A*z;  Pt =@(z) A'*z;
elseif      isobject(A)         P =@(z) A*z;  Pt =@(z) A'*z;
elseif      isa(A,'function_handle') 
    try
        if          isa(Pt,'function_handle'); P=A;
        else        error('If P is a function handle, Pt also needs to be a function handle. Exiting.'); end
    catch error('If P is a function handle, Pt needs to be specified. Exiting.'); end
else        error('P is of unsupported type. Use matrix, function_handle or object. Exiting.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Do we start from zero or not?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if initial_given ==1;
    if sum(s_initial(s_initial~=0)<lambda)
        display('Initial vector has non-zero values smaller than lambda. Setting these to zero.')
    end
    s                   =   s_initial;
    s(abs(s)<lambda)    =   0;
    Residual            =   x-P(s);
    oldERR      = Residual'*Residual/n;
else
    s_initial   = zeros(m,1);
    Residual=x;
    s=s_initial;
    oldERR      = sigsize;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Random Check to see if dictionary norm is below 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        x_test=randn(m,1);
        x_test=x_test/norm(x_test);
        nP=norm(P(x_test));
        if abs(mu*nP)>1;
            display('WARNING! Algorithm likely to become unstable.')
            display('Use smaller step-size or || P ||_2 < 1.')
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Main algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
   display('Main iterations...') 
end
tic
t=0;
done = 0;
iter=1;

while ~done

    s                   =   s + mu * Pt(Residual);
    s(abs(s)<lambda)    =   0;
    Residual            =   x-P(s);
        
    ERR=Residual'*Residual/n;
     if comp_err
         err_mse(iter)=ERR;
     end
     
     if comp_time
         iter_time(iter)=toc;
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Are we done yet?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
         if comp_err && iter >=2
             if ((err_mse(iter-1)-err_mse(iter))/sigsize<STOPTOL);
                done = 1; 
             elseif verbose && toc-t>10
                display(sprintf('Iteration %i. --- %i mse change',iter ,(err_mse(iter-1)-err_mse(iter))/sigsize)) 
                t=toc;
             end
         else
             if ((oldERR - ERR)/sigsize < STOPTOL);
                done = 1; 
             elseif verbose && toc-t>10
                display(sprintf('Iteration %i. --- %i mse change',iter ,(oldERR - ERR)/sigsize)) 
                t=toc;
             end
         end
         
     % Also stop if residual gets too small or maxIter reached
     if comp_err
         if err_mse(iter)<1e-16
             display('Stopping. Exact signal representation found!')
             done=1;
         end
     elseif iter>1 
         if ERR<1e-16
             display('Stopping. Exact signal representation found!')
             done=1;
         end
     end

     if iter >= MAXITER
         display('Stopping. Maximum number of iterations reached!')
         done = 1; 
     end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    If not done, take another round
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     if ~done
        iter=iter+1; 
        oldERR=ERR;
     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Only return as many elements as iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout >=2
    err_mse = err_mse(1:iter);
end
if nargout ==3
    iter_time = iter_time(1:iter);
end
if verbose
   display('Done') 
end
