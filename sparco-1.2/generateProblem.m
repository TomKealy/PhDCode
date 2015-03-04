function data = generateProblem(varargin)
%GENERATEPROBLEM  The main interface to individual problems.
%
%   GENERATEPROBLEM(K) creates problem structure for problem K.
%
%   GENERATEPROBLEM('list') returns a list of valid problem numbers.
%
%MATLAB SPARCO Toolbox.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: generateProblem.m 1040 2008-06-26 20:29:02Z ewout78 $

% Parse default arguments
[flags,varargin]= parseOptions(varargin,{'list','update-all','help','lookup','version'},{});
opts = parseDefaultOpts(varargin{:});

% Return version
if flags.version
   data = 'Sparco version 1.1.2'; return;
end

% Extract problem index
data = []; index = []; errstr = [];
if ~isempty(varargin) && (isscalar(varargin{1}) || ischar(varargin{1}))
   index = varargin{1}; varargin = varargin(2:end);
   
   % Look up problem index from name
   if ischar(index)
      index = lower(index);
      idx = generateProblem('list');
      for i=1:length(idx)
         % Manually construct problem call for performance
         name = eval(sprintf('prob%03d(''getname'');',idx(i)));
         if strcmp(index,lower(name)), index = idx(i); break; end;
      end
      
      if flags.lookup
        if ~isscalar(index), index = []; end;
        data = index; return;
      end
      
      if ~isscalar(index)
         disp(sprintf('Problem `%s'' not found',index));
         return
      end
   end   
end


% Handle the list request
if flags.list
   files = dir([opts.problempath 'prob*.m']);
   for i=1:length(files)
       idx = sscanf(files(i).name,'prob%d.m');
       if isnumeric(idx), data(end+1) = idx; end;
   end
   return;
end

% Handle update-all request
if flags.update_all
   problems = generateProblem('list',varargin{:});
   for i=1:length(problems)
       fprintf('Updating problem %3d . . . ', problems(i));
       data = generateProblem(problems(i),'update',varargin{:});
       fprintf('Done\n');
   end
   data = problems;
   return;
end

% Check index value
if isempty(index) && flags.help
   help generateProblem;
   return
elseif isempty(index) || (index <= 0)
   errstr = 'Problem index must be a positive scalar or valid name';
elseif ~exist(sprintf('%sprob%03d.m',opts.problempath,index),'file');
   errstr = sprintf(['Problem %d does not exists, use generateProblem', ...
                     '(''list'')\nto see all valid problems.\n'], index);
end

% Return error if needed
if ~isempty(errstr) && ~flags.lookup
   error(errstr);
   return;
elseif flags.lookup
   if ~isempty(errstr), index = []; end;
   data = index; return;
end


% Give help information on the problem
if flags.help
  help(sprintf('prob%03d',index)) 
  return
end

% Generate the problem and set info.sparcoID
[data] = eval(sprintf('prob%03d(varargin{:});',index));
if isfield(data,'info'), data.info.sparcoID = index; end;