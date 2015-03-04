function x = solveTV(M,B,TV,b,x0,options)
% SOLVETV   Demo total-variation solver
%
%  Function [X] = SOLVETV(M,B,TV,b,X0,OPTIONS), solves
%
%  Minimize ||MBx - b||^2 + gammaL1 * ||W_L1 * x     ||_p
%                         + gammaTV * ||W_TV * TV(Bx)||_q
%
%  where M is the measurement matrix and B is the sparsity basis,
%  TV is the total variation operator and W1 and W2 are diagonal
%  weight matrices for the sparsity and Total Variation terms
%  respectively.
%
%  To solve the above problem a nonlinear conjugate gradient method
%  as described in the references is used. The implementation of
%  this method in SparseMRI by M. Lustig was used as an example.
%
%  References
%     J. Nocedal and S.J. Wright, Numerical Optimization, Springer
%        series in operations research, first  edition, 1999.
%        (pp. 120-121)
%
%     M. Lustig, D.L. Donoho, and J.M. Pauly, Sparse MRI: The
%         application of compressed sensing for rapid MR imaging,
%         Submitted to Magnetic Resonance in Medicine, 2007.
%
%     SparseMRI: http://www.stanford.edu/~mlustig/SparseMRI.html

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: solveTV.m 1040 2008-06-26 20:29:02Z ewout78 $


%----------------------------------------------------------------------
% Grab input options and set defaults where needed. 
%----------------------------------------------------------------------
gradTol  = options.gradTol;
maxIts   = options.maxIts;
maxLSIts = options.maxLSIts;
p        = options.pNorm;
q        = options.qNorm;
alpha    = options.alpha;
beta     = options.beta;
mu       = options.mu;
gammaLp  = 1;
gammaTV  = 1;
weightLp = (options.weightLp(:)');
weightTV = (options.weightTV(:)' );


%----------------------------------------------------------------------
% Initialize local variables.
%----------------------------------------------------------------------
iter   = 0;
x      = x0;
stat   = 0;
fid    = 0;

% Total-variation and p-Norm terms and flags
fLp    = 0;
gLp    = 0;
fTV    = 0;
gTV    = 0;
flagLp = ((gammaLp ~= 0) && (any(weightLp ~= 0)));
flagTV = ((gammaTV ~= 0) && (any(weightTV ~= 0)));

% Exit conditions (constants).
EXIT_OPTIMAL       = 1;
EXIT_ITERATIONS    = 2;

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
logB = ' %5i  %13.7e  %13.7e';
logH = ' %5s  %13s  %13s\n';
disp(sprintf(logH,'Iter','Objective','gNorm'));


% Compute gradient and objective information
[f,g] = computeInfo(x);
dx = -g;
gNorm = sqrt(g(:)'*g(:));


%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------
while 1
   
   %-------------------------------------------------------------------
   % Test exit conditions.
   %-------------------------------------------------------------------
   if (iter  >= maxIts), stat = EXIT_ITERATIONS; end;
   if (gNorm < gradTol), stat = EXIT_OPTIMAL;    end;

   %-------------------------------------------------------------------
   % Print log and act on exit conditions.
   disp(sprintf(logB,iter,f,gNorm));
   
   if stat, break; end % Act on exit conditions.

   %===================================================================
   % Iterations begin here.
   %===================================================================
   iter = iter + 1;

   % Backtracking line-search
   [xNew,fNew,gNew,lnErr] = linesearch(x,f,g,dx,alpha,beta,maxLSIts);
   if lnErr
      stat = EXIT_LINE_ERROR;
      break;
   end
   
   % Compute gradient norm
   gNormNew = sqrt(gNew' * gNew);
    
   % Compute gamma
   gamma = (gNormNew / (gNorm+eps)).^2;

   % Update x, dx, f, g and gNorm
   dx    = -gNew + gamma * dx;
   x     = xNew;
   f     = fNew;
   g     = gNew;
   gNorm = gNormNew;
end

% Print final output.
switch (stat)
   case EXIT_OPTIMAL
      disp(sprintf('\n EXIT -- Optimal solution found\n'));
   case EXIT_ITERATIONS
      disp(sprintf('\n ERROR EXIT -- Too many iterations\n'));
   case EXIT_LINE_ERROR
      disp(sprintf('\n ERROR EXIT -- Linesearch error (%i)\n',lnErr));
   otherwise
      error('Unknown termination condition\n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f,g] = computeInfo(x)
   z    = B(x,1);
   Mz   = M(z,1);
   Tz   = TV(z,1);
   Tzw  = Tz.*weightTV;
   xw   = x.*weightLp;
   Tz2w = Tzw.*conj(Tzw);
   x2w  = xw.*conj(xw);

   % Compute the objective
   Mzb = Mz - b;
   fRes = Mzb' * Mzb;
   if flagLp, fLp = sum(power(x2w + mu, p/2));  end
   if flagTV, fTV = sum(power(Tz2w(:) + mu, q/2));  end
   f = fRes + gammaLp * fLp + gammaTV * fTV;

   % Compute the gradient
   gRes = 2 * B(M(Mzb,2),2);
   if flagLp,
     gLp = p *     (weightLp.*(xw.*power(x2w   + mu,p/2-1)));
   end
   if flagTV,
     gTV = q * B(TV(weightTV.*(Tzw.*power(Tz2w + mu,q/2-1)),2),2);
   end;
   g = gRes + gammaLp * gLp + gammaTV * gTV;
end


function [xNew,fNew,gNew,err] = linesearch(x,f,g,dx,alpha,beta,maxIts)
   step = 1; err = 1;
   while (maxIts > 0)
     xNew = x + step * dx;
     [fNew,gNew] = computeInfo(xNew);

     % Check Armijo condition
     if (fNew <= f + alpha*step*real(g'*dx));
       err = 0; break;
     end;

     % Reduce steplength and number of iterations left
     step   = step * beta;
     maxIts = maxIts - 1;
   end
end


end % solveTV
