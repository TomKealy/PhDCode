function exTV
% Illustration of Total Variation minimization
%
% References
%     M. Lustig, D.L. Donoho, and J.M. Pauly, Sparse MRI: The
%         application of compressed sensing for rapid MR imaging,
%         Submitted to Magnetic Resonance in Medicine, 2007.
%
%     SparseMRI: http://www.stanford.edu/~mlustig/SparseMRI.html

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: exTV.m 1040 2008-06-26 20:29:02Z ewout78 $

  % Create the data
  data = generateProblem(502);

  % Create total difference 'operator'
  TV = opDifference(data.signalSize);

  % Set solver parameters
  opts.maxIts           = 10;
  opts.maxLSIts         = 150;
  opts.gradTol          = 1e-30;

  opts.weightTV         = 0.001;
  opts.weightLp         = 0.01;
  opts.pNorm            = 1;
  opts.qNorm            = 1;
  opts.alpha            = 0.01;
  opts.beta             = 0.6;
  opts.mu               = 1e-12;

  % Give instructions to the user.
  fprintf('This script calls "solveTV" five times,\n');
  fprintf('each with a maximum of 10 iterations.\n');
  fprintf('Ignore the messages "ERROR EXIT"\n');
  input('Press "Return" to continue.');
  
  % Solve
  x = randn(prod(data.signalSize),1);
  for i=1:5
    x = solveTV(data.M, data.B, TV, data.b, x, opts);
    y = data.reconstruct(x);
    imagesc(abs(y));
    title(sprintf('Iteration %d',i));
    drawnow;
  end
  
end % function exTV

