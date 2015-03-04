function data = completeOps(data)

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: completeOps.m 1040 2008-06-26 20:29:02Z ewout78 $

operators = {};
flagM = 0; if isfield(data,'M'), flagM = 1; end;
flagB = 0; if isfield(data,'B'), flagB = 1; end;

if (~flagM) && (~flagB)
   error('At least one of the operators M or B has be to given.');
end

% Define measurement matrix if needed
if ~flagM
  info   = data.B([],0);
  data.M = opDirac(info{1});
else
  operators{end+1} = data.M;
end

% Define sparsity basis if needed
if ~flagB
  info   = data.M([],0);
  data.B = opDirac(info{2});
else
  operators{end+1} = data.B;
end

% Define operator A if needed
if ~isfield(data,'A')
  if (length(operators) > 1)
     data.A = opFoG(operators{:});
  else
     data.A = operators{1};
  end
end

% Define empty solution if needed
if ~isfield(data,'x0')
  data.x0 = [];
end

% Define the operator size and string
opInfo       = data.A([],0);
data.sizeA   = [opInfo{1},opInfo{2}];
opInfo       = data.B([],0);
data.sizeB   = [opInfo{1},opInfo{2}];
opInfo       = data.M([],0);
data.sizeM   = [opInfo{1},opInfo{2}];
data.op.strA = opToString(data.A);
data.op.strB = opToString(data.B);
data.op.strM = opToString(data.M);

% Get the size of the desired signal
if ~isfield(data,'signalSize')
   if ~isfield(data,'signal')
     error(['At least one of the fields signal ', ...
            'or signalSize must be given.']);
   end
   data.signalSize = size(data.signal);
end

% Reconstruct signal from sparse coefficients
if ~isfield(data,'reconstruct')
   data.reconstruct = @(x) reshape(data.B(x,1),data.signalSize);
end

% Reorder the fields (sort alphabetically)
m = fieldnames(data);
m = sort(m);
dataReorder = struct();
for i=1:length(m)
  eval(sprintf('dataReorder.%s = data.%s;',m{i},m{i}));
end

data = dataReorder;
