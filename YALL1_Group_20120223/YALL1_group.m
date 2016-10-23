function [x,Out] = YALL1_group(A,b,groups,varargin)
%
%  This is a solver for
%
%  (1) group-sparse basis pursuit model:
%
%      minimize     sum_i w_i*||x_{gi}||_2
%      subject to   Ax = b,
%                   x >= 0 (optional),
%
%      where 
%      - w_i is the weight for the i-th group; 
%      - gi is the index set of the i-th group;
%      - groups may overlap.
%
%  (2) jointly-sparse basis pursuit model:
%
%      minimize     sum_i w_i*||X(i,:)||_2
%      subject to   AX = B [or A{l}X(:,l) = B(:,l), for l = 1,2,...,L]
%                   X >= 0 (optional),
%
%      where
%      - w_i is the weight for the i-th row;
%      - X(i,:) is the i-th row of matrix X;
%      - the sensing matrix can either be the same A for each 
%        column of X, or be different A{l}, l= 1,2,...L.
% 
% ----------------------------------------------------------
%  Authors: Wei Deng, Wotao Yin, Yin Zhang                  
%           Rice University                                                                                            
% ----------------------------------------------------------
%
% ===================== Required inputs =====================
%
%  A -- multiple types of A can be accepted:
%       1) an m x n matrix, 
%       2) a cell array of m x n matrices (for joint-sparse model with
%          different sensing matrices)
%       3) a structure with fields:
%          a) A.times: a function handle for A*x
%          b) A.trans: a function handle for A'*x
%          c) A.invIpAAt: a function handle for (beta1*eye(m)+beta2*(A*A'))\x,
%             which is required when
%             - primal solver is to be used, AND
%             - A is non-orthonormal, AND
%             - exact linear system solving is to be performed.
%
%          d) A.invAAt: a function handle for (A*A')\x, which is required when
%             - dual solver is to be used, AND
%             - A is non-orthonormal, AND
%             - exact linear system solving is to be performed.
%
%  b -- an m-entry vector for the group-sparse model, or
%       an m x l matrix for the jointly-sparse model
%
%  groups -- 1) For non-overlapping group-sparse model: an n-entry vector  
%               whose i-th entry is the group number of x_i; 
%            2) For overlapping group-sparse model: a cell array whose   
%               i-th cell constains the indices of i-th group in a column 
%               vector;               
%            3) For jointly-sparse model: [].
%
% ===================== Optional inputs =====================
%
%  'StopTolerance' -- stopping tolerance
%
%  'GrpWeights' -- weights for the groups
% 
%  'overlap' -- true for overlapping groups
%
%  'Nonnegative' -- true for imposing nonnegativity constraints
%
%  'nonorthA' -- true for A with non-orthonormal rows; false otherwise
%
%  'ExactLinSolve' -- true for solving linear systems exactly;
%                     false for taking a gradient descent step
%
%  'QuadPenaltyPar' -- nonnegative penalty parameters
% 
%  'StepLength' -- nonnegative step lengths
%
%  'maxIter' -- maximum number of iterations
%
%  'xInitial' -- an initial estimate of the solution
%
%  'Solver' -- 1 for primal-based solver; 
%              2 for dual-based solver
%
%  'Continuation' -- false for off; true for performing continuation on the
%                    penalty parameters: it allows small initial penalty 
%                    parameters for constraint violations, which lead to 
%                    faster initial convergence, and it increases those 
%                    parameters whenever the violation reduction slows down.
%                    It leads to overall speedups in most cases.
%                    
%                    We use the following continuation scheme: 
%                    multiply the penalty parameters by a factor c (c > 1) 
%                    if ||R||_2 > alpha*||R_prev||_2, where 0 < alpha < 1 
%                    is a parameter, R and R_prev denote the constraint 
%                    violations at the current and previous iterations, 
%                    respectively.
%
%  'ContParameter' -- the parameter alpha (0 < alpha < 1) in the 
%                     continuation scheme
%
%  'ContFactor' -- the factor c (c > 1) in the continuation scheme
%
% ===================== Outputs =====================
%
%  x -- last iterate (hopefully an approximate solution)
%
%  Out -- a structure with fields:
%         Out.status  -- exit information
%         Out.iter    -- # of iterations taken
%         Out.cputime -- solver CPU time
%
% ----------------------------------------------------------

% Test for number of required parametres
if (nargin-length(varargin)) ~= 3
    error('Wrong number of required parameters');
end

m = size(b,1); % # of measurements

%--------------------------------------------------------------
% Set parameters to their defaults  
%--------------------------------------------------------------
opts.tol = 1e-6;
opts.nonneg = false;
opts.exact = false;
opts.maxit = 10*m;
opts.solver = 2;
opts.continuation = false;
opts.contpar = 0.9;
opts.contfactor = 1.2;

%--------------------------------------------------------------
% Set parameters to user specified values
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Options should be given in pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STOPTOLERANCE'
                opts.tol = varargin{i+1};
            case 'GRPWEIGHTS'
                opts.weights = varargin{i+1};
            case 'OVERLAP'
                opts.overlap = varargin{i+1};
            case 'NONNEGATIVE'
                opts.nonneg = varargin{i+1};
            case 'NONORTHA'
                opts.nonorth = varargin{i+1};
            case 'EXACTLINSOLVE'
                opts.exact = varargin{i+1};
            case 'QUADPENALTYPAR'
                opts.beta = varargin{i+1};
            case 'STEPLENGTH'                
                opts.gamma = varargin{i+1};
            case 'MAXITER'
                opts.maxit = varargin{i+1};
            case 'XINITIAL'
                opts.xInit = varargin{i+1};
            case 'SOLVER'
                opts.solver = varargin{i+1};
            case 'CONTINUATION'
                opts.continuation = varargin{i+1};
            case 'CONTPARAMETER'
                opts.contpar = varargin{i+1};
            case 'CONTFACTOR'
                opts.contfactor = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

%--------------------------------------------------------------
% Set the default ADM parameters if not specified
%--------------------------------------------------------------
if opts.solver == 1  % primal solver
    if ~isfield(opts,'beta')
        opts.beta = [1,1]/mean(abs(b(:)));
    elseif length(opts.beta) ~= 2
        error('Option ''QuadPenaltyPar'' must be a 2-vector for the primal solver.');       
    end
    if ~isfield(opts,'gamma')
        opts.gamma = [1.618, 1.618];
    elseif length(opts.gamma) ~= 2
        error('Option ''StepLength'' must be a 2-vector for the primal solver.');
    end
elseif opts.solver == 2  % dual solver
    if ~isfield(opts,'beta')
        opts.beta = 1*mean(abs(b(:)));
    elseif length(opts.beta) ~= 1
        error('Option ''QuadPenaltyPar'' must be a scalar for the dual solver.');
    end
    if ~isfield(opts,'gamma')
        opts.gamma = 1.618;
    elseif length(opts.gamma) ~= 1
        error('Option ''StepLength'' must be a scalar for the dual solver.');
    end
else
    error('Option ''Solver'' must be 1 or 2')
end

if ~isfield(opts,'overlap')
    if iscell(groups)
        opts.overlap = true;
    else
        opts.overlap = false;
    end
end

if opts.overlap, 
    opts.nonorth = true; opts.exact = false;
end

% Define linear operators
[A,At,R] = linear_operators(A,b,opts);

%--------------------------------------------------------------
% Read joint/group sparsity information
%--------------------------------------------------------------
if isempty(groups)  % jointly-sparse model
    solver = opts.solver+2;
    nG = length(At(zeros(m,1)));  % # of rows of the solution matrix
elseif opts.overlap  % overlapping group-sparse model
    solver = opts.solver+4;
    nG = length(groups); % # of groups
    n = length(At(zeros(m,1)));
    % Preprocess the groups
    gcat = cat(1,groups{:});
    N = length(gcat);
    gidx = unique(gcat);
    if min(gidx)~=1 || max(gidx)~=n
        error('Please normalize the group numbering');
    end
    gmat = sparse(nG,N);
    ii = 1;
    for i = 1:nG
        gmat(i,ii:ii+length(groups{i})-1) = 1;
        ii = ii + length(groups{i});
    end
    for i = 1:nG
        groups{i} = i*ones(length(groups{i}),1);
    end
    groups = cat(1,groups{:});   
    G = sparse(N,n);
    for i = 1:N
        G(i,gcat(i)) = 1;
    end  
else  % non-overlapping group-sparse model
    solver = opts.solver;
    n = length(groups);
    gidx = unique(groups);
    nG = length(gidx);  % # of groups   
    if min(gidx)~=1 || max(gidx)~=nG
        error('Please normalize the group numbering');
    end    
    % Preprocess the groups
    gmat = sparse(nG,n);
    for i=1:nG
        gmat(i,groups == i) = 1;
    end
end

if ~isfield(opts,'weights')
    opts.weights = ones(nG,1);  % default weights
end

% Check whether the rows of A are orthonormal if not specified
if ~isfield(opts,'nonorth')
    opts.nonorth = check_orth(A,At,b);
end

% If A has orthonormal rows, always solve the subproblems exactly
if ~opts.nonorth, opts.exact = true; end

tic;  % start the clock
%--------------------------------------------------------------
% Call the ADM solver
%--------------------------------------------------------------
addpath([fileparts(mfilename('fullpath')) '/main_solvers']);
addpath([fileparts(mfilename('fullpath')) '/subprob_solvers']);
switch solver  
    case 1  % primal solver for the non-overlapping group-sparse model
        [x,Out] = group_primal_solver(A,At,R,b,groups,gmat,opts);
    case 2  % dual solver for the non-overlapping group-sparse model
        [x,Out] = group_dual_solver(A,At,R,b,groups,gmat,opts);
    case 3  % primal solver for the jointly-sparse model
        [x,Out] = joint_primal_solver(A,At,R,b,opts);
    case 4  % dual solver for the jointly-sparse model
        [x,Out] = joint_dual_solver(A,At,R,b,opts);
    case 5  % primal solver for the overlapping group-sparse model
        [x,Out] = overlap_primal_solver(A,At,b,groups,gmat,G,opts);
end
Out.cputime = toc;

end

%% ===============================================================
function [A,At,R] = linear_operators(A0,b,opts)
%-----------------------------------------------------------------
%  Define linear operators:
%
%  A  -- a function handle for A*x
%  At -- a function handle for A'*x
%  R (optional) -- a function handle for 
%                  1) (beta1*eye(m)+beta2*(A*A'))\x for primal solver
%                  2) (A*A')\x for dual solver
%
%------------------------------------------------------------------
R = [];
if isnumeric(A0)
    if opts.nonorth && opts.exact
        if opts.solver == 1
            m = size(b,1);
            R = inv(opts.beta(1)*speye(m)+opts.beta(2)*(A0*A0'));
        else
            R = inv(A0*A0');
        end
        R = @(x) R*x;
    end
    A  = @(x) A0*x;
    At = @(x) (x'*A0)';
elseif isstruct(A0) && isfield(A0,'times') && isfield(A0,'trans');
    A  = A0.times;
    At = A0.trans;
    if opts.nonorth && opts.exact
        if opts.solver == 1 && isfield(A0,'invIpAAt')
            R = @(x) A0.invIpAAt(x,opts.beta);
        elseif opts.solver == 2 && isfield(A0,'invAAt')
            R = @(x) A0.invAAt(x);
        else
            error('Field .invIpAAt or .invAAt is missing.')
        end
    end
elseif iscell(A0) && isnumeric(A0{1})
    if opts.nonorth && opts.exact
        Rcell = cell(size(A0));
        if opts.solver == 1
            m = size(b,1);
            for l = 1:length(A0)	% parfor ready
                Rcell{l} = inv(opts.beta(1)*speye(m)+opts.beta(2)*(A0{l}*A0{l}'));
            end
        else
            for l = 1:length(A0)	% parfor ready
                Rcell{l} = inv(A0{l}*A0{l}');
            end
        end
        R = @(x) cell_m_matrix(Rcell,x);
    end
    A  = @(x) cell_m_matrix(A0,x);
    At = @(x) cellT_m_matrix(A0,x);
else
    error('A must be a matrix or a structure with fields .times and .trans');
end
end

% function for defining A*X for cell array A
function res = cell_m_matrix(AA,XX)
    % compute res = [ ... AA{l}*XX(:,l) ... ]
    MM = size(AA{1},1);
    LL = size(XX,2);
    res = zeros(MM,LL);
    for ll = 1:LL  % parfor ready
        res(:,ll) = AA{ll}*XX(:,ll);
    end
end
% function for defining A'*Y for cell array A
function res = cellT_m_matrix(AA,YY)
    % compute res = [ ... AA'{l}*YY(:,l) ... ]
    NN = size(AA{1},2);
    LL = size(YY,2);
    res = zeros(NN,LL);
    for ll = 1:LL  % parfor ready
        res(:,ll) = (YY(:,ll)'*AA{ll})';
    end
end


%% ===============================================================
function nonorth = check_orth(A,At,b)
%------------------------------------------------
% check whether the rows of A are orthonormal
%------------------------------------------------
nonorth = false;
s1 = randn(size(b,1));
s2 = A(At(s1));
err = norm(s1-s2)/norm(s1);
if err > 1.e-12; nonorth = true; end
end