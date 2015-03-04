function checkProblems(list)
%CHECKPROBLEMS  Generate problems and check operator A
%
%   CHECKPROBLEMS(LIST) creates the default instances of all
%   problems given in LIST. For each problem, operator A is checked
%   by the DOTTEST program to detect mistakes.
%
%   See also DOTTEST.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: checkProblems.m 1040 2008-06-26 20:29:02Z ewout78 $

if (nargin < 1)
  list = generateProblem('list');
end

cws = get(0,'CommandWindowSize');

for i = list
   % Generate the test problem
   try
      P = generateProblem(i);
      status = 0;
   catch
      P.info.title = '';
      status = -1;
   end   
   
   width = max(0,cws(1)-length(P.info.title)-30);
   dots  = repmat('. ', 1,floor(width/2));
   fprintf('Checking problem %3d: %s %*s', i, P.info.title, width, ...
           dots);
   
   % Check if problem generation was ok
   if status == -1
     err = lasterror; msg = err.message;
     fprintf('Failed\n');
     fprintf('%22sWARNING: Failed to instantiate problem %d\n','', i);
     fprintf('%22s         %s\n', '', msg);
     continue;
   end
   
   % Perform dot test on operator P.A
   try
      status = dottest(P.A,5,'quiet');
   catch
      status = -1;
   end
   
   if status == 1
     fprintf('Failed\n');
     fprintf('%22sWARNING: Dot test failed on operator A!\n', '', i);
   elseif status == -1
     err = lasterror; msg = err.message;
     fprintf('Failed\n');
     fprintf('%22sWARNING: %s\n','',strrep(msg,char(10),': '));
   else
     fprintf('OK\n');
   end
end
