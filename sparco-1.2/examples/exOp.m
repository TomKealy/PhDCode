function exOp
% Example illustrating the use of operators.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: exOp.m 1040 2008-06-26 20:29:02Z ewout78 $

  % Set problem size
  k = 20;  % Number of non-zero coefficients
  n = 128; % Coefficient length
  m = 80;  % Number of measurements

 
  % Set up sparse coefficients
  p = randperm(n);
  x0= zeros(n,1); x0(p(1:k)) = randn(k,1);

  % Set up operators
  disp(['Setting up partial DCT measurement ' ...
        'and Haar sparsity operators']);
  p = randperm(n);
  M = opFoG(opRestriction(n,p(1:m)),opDCT(n));
  B = opHaar(n);
  A = opFoG(M,B);
  
  % Create operator class for B
  disp('Turning Haar operator into a class');
  Haar = classOp(B);  

  % Compute observed signal b
  disp('Computing b = M*B*x0');
  b = M(Haar * x0,1);
  
  % Solve the problem
  input('Press <RETURN> to solve the problem . . .');
  [xs,r,g,info] = spgl1(A,b,0,0,[]);
  
  figure(1);
  plot(1:n,xs,'b',1:n,x0,'ro');
  legend('Solution found','Original coefficients');
  
  figure(2);
  y = Haar * xs; % Signal from sparse Haar coefficients
  plot(1:n,y,'b-',1:n,Haar * x0,'ro');
  legend('Reconstructed signal','Original signal');
  
  % Quick check on operator
  disp('Corroborate correctness of operator A by checking that')
  disp('<Ax,y> = <x,A* y> for 100 random vectors x and y.');
  dottest(A,100);
  
  % Extract matrix entries and plot
  disp('Extracting matrix from operator (plot 3)');
  M = opToMatrix(A);
  figure(3); imagesc(M);

end % function exOp
