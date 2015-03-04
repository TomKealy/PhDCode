function status = dottest(op,n,mode)
%DOTTEST  Detect errors in the implementation of an operator
%
%   DOTTEST(OP,N,MODE) generates N random real and complex vectors
%   X and Y, and verifies that OP(X,1)'*Y = X'*OP(Y,2) within a
%   tolerance of 1E-10. This can help detect errors in the
%   operator; it canot be used to guarantee correctness. The
%   function returns 0 when the test succeeded and -1 if it
%   failed. Depending on the domain of the operator (i.e., the
%   complexity flags obtained from OP([],0)), real or complex 
%   vector checks may be omitted. MODE can be set to 'quiet' to
%   suppress output. Parameters N and MODE are optional and are set
%   to 100 and '' respectively by default.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: dottest.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 2), n = 100; end;
if (nargin < 3), mode = ''; end;

tol    = 1e-10;
info   = op([],0);
status = 0;
quiet  = 0;

if strcmp(mode,'quiet'), quiet = 1; end;

% Print testing info
if ~quiet
   name = opToString(op);
   fprintf('Performing dot test on operator: %s\n', name);
end

% See if we want to test real vectors
domain = info{3};
if (domain(1)==0 && domain(3)==0)
  testReal = 1;
else
  testReal = 0;
end

if (domain(2)==1) && (domain(4)==1)
  testComplex = 1;
else
  testComplex = 0;
end

if testComplex
  complexPassed = 0; complexErrmax = 0; complexRatio = [inf,0];
  for i=1:n
    x  = randn(info{2},1) + sqrt(-1)*randn(info{2},1);
    y  = randn(info{1},1) + sqrt(-1)*randn(info{1},1);
    z1 = op(x,1)' * y;
    z2 = x' * op(y,2);
    err = abs(z1 - z2);
    if (err > complexErrmax), complexErrmax = err; end;
    if (err < tol),
      complexPassed = complexPassed + 1;
    else
      ratio = abs(z1) / abs(z2);
      if ratio < complexRatio(1), complexRatio(1) = ratio; end;
      if ratio > complexRatio(2), complexRatio(2) = ratio; end;
      status = -1;
    end
  end
end

if testReal
  realPassed = 0; realErrmax = 0; realRatio = [inf,0];
  for i=1:n
    x = randn(info{2},1);
    y = randn(info{1},1);
    z1 = op(x,1)' * y;
    z2 = x' * op(y,2);
    err = abs(z1 - z2);
    if (err > realErrmax), realErrmax = err; end;
    if (err < tol),
      realPassed = realPassed + 1;
    else
      ratio = abs(z1) / abs(z2);
      if ratio < realRatio(1), realRatio(1) = ratio; end;
      if ratio > realRatio(2), realRatio(2) = ratio; end;
      status = -1;
    end
  end
end


if ~quiet
  if testComplex
     if complexPassed < n
        fprintf('Complex: FAILED on %d out of %d tests\n', ...
                n-complexPassed, n);
        fprintf('%8s maximum absolute difference of %13.9e\n', ...
                '',complexErrmax);
        fprintf('%8s ratio between %13.9e and %13.9e\n', ...
                '',complexRatio(1),complexRatio(2));
     else
        fprintf('Complex: PASSED!\n');
     end
   else
     fprintf('Complex: -\n');
   end

   if testReal
     if realPassed < n
        fprintf('Real  : FAILED on %d out of %d tests\n', ...
                n-realPassed, n);
        fprintf('%8s maximum absolute difference of %13.9e\n', ...
                '',realErrmax);
        fprintf('%8s ratio between %13.9e and %13.9e\n', ...
                '',realRatio(1),realRatio(2));
     else
        fprintf('Real   : PASSED!\n');
     end
   else
     fprintf('Real   : -\n');
   end
end
