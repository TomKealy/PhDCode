function [sol, numIters, activationHist, v] = SolvePFP(A, y, N, algType, maxIters, verbose, OptTol)
% SolvePFP: Implements the Polytope Faces Pursuit algorithm
% Usage
%	[sols, numIters, activationHist] = SolvePFP(A, y, N, algType, maxIters,
%	verbose, OptTol)
% Input
%	A           An nxN matrix, with rank(A) = min(N,n).
%	y           vector of length n.
%   N           length of solution vector. 
%   algType     'lars' for the Lars algorithm, 
%               'pfp' for Polytope Faces Pursuit (default).
%               Add prefix 'nn' (i.e. 'nnlars' or 'nnpfp') to add a
%               non-negativity constraint (omitted by default)
%	maxIters    maximum number of Lars iterations to perform. If not
%               specified, runs to stopping condition (default)
%   solFreq     if =0 returns only the final solution, if >0, returns an 
%               array of solutions, one every solFreq iterations (default 0). 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default)
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	sols           solution of the PFP algorithm
%	numIters       Total number of steps taken
%   activationHist Array of indices showing elements entering and 
%                  leaving the solution set
% Description
%   SolvePFP implements the Polytope Faces Pursuit algorithm, to 
%   solve the problem 
%      min || x ||_1  s.t. Ax = y
%   via its dual
%      max y'*v s.t. A'*v <= 1
%   The implementation implicitly factors the active set matrix A(:,I)
%   using Choleskly updates. 
% References
%   M.D. Plumbley, "Recovery of Sparse Representations by Polytope Faces
%   Pursuit", In Proceedings of the 6th International Conference on
%   Independent Component Analysis and Blind Source Separation (ICA 2006), 
%   Charleston, SC, USA, 5-8 March 2006, LNCS 3889, pp 206-213. 
%   Springer-Verlag, Berlin, 2006
% See Also
%   SolveOMP, SolveBP, SolveStOMP, SolveLasso
%

n = length(y);

if nargin < 8,
    OptTol = 1e-5;
end
if nargin < 7,
    verbose = 0;
end
if nargin < 6,
    solFreq = 0;
end
if nargin < 5,
    maxIters = 10*n;
end
if nargin < 4,
    algType = 'pfp';
end
if nargin < 3,
    N = size(A,2);
end

switch lower(algType)
    case 'lars'
        isLars = 1; nonNegative = 0;
    case 'pfp'
        isLars = 0; nonNegative = 0;
    case 'nnlars'
        isLars = 1; nonNegative = 1;
    case 'nnpfp'
        isLars = 0; nonNegative = 1;
end

if ~nonNegative
    A = [A -A];
    N = 2*N;
end

% Global variables for linsolve function
global opts opts_tr errTol
opts.UT = true; 
opts_tr.UT = true; opts_tr.TRANSA = true;
errTol = 1e-5;

k = 0;
v = zeros(n,1);
res = y;
R_I = [];
activeSet = [];
activationHist = activeSet;
collinearIndices = [];
done = 0;
while  ~done
    k = k+1;
    
    % Find next face
    inactiveSet = 1:N;
    inactiveSet(activeSet) = 0;
    inactiveSet(collinearIndices) = 0;
    inactiveSet = find(inactiveSet > 0);

    Ares = A(:,inactiveSet)'*res;
    Av = A(:,inactiveSet)'*v;

    if max(Ares) > OptTol
        posInd = find(Ares > 0);
        [alpha_p,face_p] = max(Ares(posInd) ./ (1 - Av(posInd)));
        newFace = inactiveSet(posInd(face_p));
    else
        k = k-1;
        done = 1;
        break;
    end

    % Add new face to active set
    % Update the Cholesky factorization of A_I
    [R_I, flag] = updateChol(R_I, n, N, A, activeSet, newFace);
    % Check for collinearity
    if (flag)
        collinearIndices = [collinearIndices newFace];
        if verbose
            fprintf('Iteration %d: Variable %d is collinear\n', iter, newIndices(j));
        end
        valid = 1;
    else
        activeSet = [activeSet newFace];
        activationHist = [activationHist newFace];
        if verbose
            fprintf('Iteration %d: Adding variable %d\n', k, newFace);
        end
        valid = 0;
    end

    % Find first active vector to violate sign constraint
    while ~valid
        % Compute the solution estimate - solve the equation (A_I'*A_I)x_I = y
        corr_I = A(:,activeSet)'*y;
        z = linsolve(R_I,corr_I,opts_tr);
        x_I = linsolve(R_I,z,opts);
        
        ind = find(x_I < 0);
        if (length(ind) > 0) && (~isLars)
            col = ind(1);
            
            if verbose
                fprintf('Iteration %d: Dropping variable %d\n', k, activeSet(ind));
            end

            % If violating element is the one just added, ignore it in the
            % next iteration to avoid cycling
            if (activeSet(ind) == newFace)
                collinearIndices = [collinearIndices ind];
            else
                collinearIndices = [];
            end

            % Downdate the Cholesky factorization of A_I
            R_I = downdateChol(R_I,col);
            activationHist = [activationHist -activeSet(ind)];
            activeSet = [activeSet(1:col-1), activeSet(col+1:length(activeSet))];
        else
            valid = 1;
        end
    end

    z = linsolve(R_I,ones(length(activeSet),1),opts_tr);
    z = linsolve(R_I,z,opts);
    v = A(:,activeSet)*z;
    y_k = A(:,activeSet)*x_I;
    res = y - y_k;

    if k >= maxIters
        done = 1;
    end
end

x = zeros(N,1);
x(activeSet) = x_I;
if ~nonNegative
    sol = x(1:N/2) - x(N/2+1:N);
    activationHist = (mod(abs(activationHist)-1,N/2)+1) .* sign(activationHist);
else
    sol = x;
end
numIters = k;

clear opts opts_tr errTol


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, flag] = updateChol(R, n, N, A, activeSet, newIndex)
% updateChol: Updates the Cholesky factor R of the matrix 
% A(:,activeSet)'*A(:,activeSet) by adding A(:,newIndex)
% If the candidate column is in the span of the existing 
% active set, R is not updated, and flag is set to 1.

global opts_tr errTol
flag = 0;

newVec = A(:,newIndex);

if length(activeSet) == 0,
    R = sqrt(sum(newVec.^2));
else
    p = linsolve(R,A(:,activeSet)'*A(:,newIndex),opts_tr);
    q = sum(newVec.^2) - sum(p.^2);
    if (q <= errTol) % Collinear vector
        flag = 1;
    else
        R = [R p; zeros(1, size(R,2)) sqrt(q)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = downdateChol(R, j)
% downdateChol: `Downdates' the cholesky factor R by removing the 
% column indexed by j.

% Remove the j-th column
R(:,j) = [];
[m,n] = size(R);

% R now has nonzeros below the diagonal in columns j through n.
% We use plane rotations to zero the 'violating nonzeros'.
for k = j:n
    p = k:k+1;
    [G,R(p,k)] = planerot(R(p,k));
    if k < n
        R(p,k+1:n) = G*R(p,k+1:n);
    end
end

% Remove last row of zeros from R
R = R(1:n,:);

