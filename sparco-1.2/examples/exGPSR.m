function exGPSR
% Example illustrating how problems generated with Sparco can be
% solved using GPSR.
%
% References:
%
%   M. Figueiredo, R. Nowak, and S. J. Wright, "Gradient projection
%      for sparse reconstruction: Application to compressed sensing
%      and other inverse  problems", To appear in IEEE Journal of
%      Selected Topics in Signal Processing: Special Issue on
%      Convex Optimization Methods for Signal Processing, 2007.
%
%   GPSR: http://www.lx.it.pt/~mtf/GPSR

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: exGPSR.m 1040 2008-06-26 20:29:02Z ewout78 $

% Generate problem 6: Piecwise cubic polynomial signal.
  P = generateProblem(6);

  A  = @(x) P.A(x,1);  % The operator
  AT = @(y) P.A(y,2);  % and its transpose.

  b  = P.b;            % The right-hand-side vector.
  
  m = P.sizeA(1);      % m is the no. of rows.
  n = P.sizeA(2);      % n is the no. of columns.
  
% Solve an L1 recovery problem:
% minimize  1/2|| Ax - b ||_2^2  +  lambda ||x||_1.
  lambda = 1000;
  x = GPSR_BB(b, A, lambda, 'AT', AT);
  
% The solution x is the reconstructed signal in the sparsity basis. 
  figure;
  plot(x);
  title('Coefficients of the reconstructed signal')
  xlabel('Coefficient')
  
% Use the function handle P.reconstruct to use the coefficients in
% x to reconstruct the original signal.
  y     = P.reconstruct(x);    % Use x to reconstruct the signal.
  yorig = P.signal;            % P.Signal is the original signal.

  figure;
  plot(1:length(y), y    ,'b', ...
       1:length(y), yorig,'r');
  legend('Reconstructed signal','Original signal');
  title('Reconstructed and original signals');
  
end % function exGPSR
