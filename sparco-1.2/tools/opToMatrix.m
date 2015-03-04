function A = opToMatrix(op,show)
%OPTOMATRIX  Generate explicit matrix form of operator
%
%   OPTOMATRIX(OP,SHOW) generates the explicit matrix form of
%   operator OP thought the multiplication of all columns of the
%   identity matrix. To see the progress, optional parameter SHOW
%   can be set to 1. This will give a progress bar with approximate
%   duration and percentage finished.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opToMatrix.m 1040 2008-06-26 20:29:02Z ewout78 $

% Get operator information
info = op([],0);
m    = info{1};
n    = info{2};

% Check show parameter
if nargin < 2 || isempty(show) || show==0
  show = 0;
else
  show = 1;
end
  

% Generate empty matrix. Besides speeding up the process this also
% ensures that the matrix fits in memory before generating the
% individual columns.
A = zeros(m,n);

tic; flagProgress = 0;
for i=1:n
  
  % Extract the i-th column of A
  ei = zeros(n,1); ei(i) = 1;
  A(:,i) = op(ei,1);
  
  % Make an estimation for total time
  if (i == 1) && show
    t = toc;
    if (t*n) > 1
      flagProgress = 1;
      h = waitbar(0,'Time left --/-- sec.');
    end;
  end
  
  % Report progress
  if flagProgress && (mod(i,10)==0)
    pct = floor(100.0 * i / double(n)); % Percentage done
    etl = ceil((n-i) * (toc / i));     % Estimated time left
    ett = ceil(n     * (toc / i));     % Estimated total time
    
    waitbar(pct/100, h, sprintf('Time left %d/%d sec.', etl, ett));
    set(h, 'Name',sprintf('%3d%%',pct));
  end
  
end

if flagProgress
  close(h)
end
