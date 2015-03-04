function op = opSum(varargin)
% OPSUM   Summate a sequence of operators into a single operator.
%
%    OPSUM(WEIGHTS,OP1,OP2,...OPn) creates an operator that
%    consists of the weighted sum of OP1, OP2, ..., OPn. An
%    alternative but equivalent way to call this function is
%    OPSUM(WEIGHTS,{OP1,OP2,...,OPn}). The WEIGHT parameter is
%    an optional N-vector.
%
%    See also opFoG, opDictionary, opStack.

%   Copyright 2008, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opSum.m 1040 2008-06-26 20:29:02Z ewout78 $

opList = varargin;

if nargin == 0
   error('At least one operator must be specified');
else
   % Get number of operators and weights
   if (nargin > 1) && isnumeric(varargin{1})
     opList = varargin(2:nargin);
     w = varargin{1}; w = w(:);
   else
     opList = varargin;
     w = [];
   end
   if iscell(opList{1}), opList = opList{1}; end;
   if isempty(w), w = ones(length(opList),1); end;
   
   % Check length of weights
   if (length(w) ~= length(opList))
      error('Weight vector must have same length as operator list');
   end
   
   % Check operator consistency and domain (real/complex)
   info = opList{1}([],0);
   m = info{1}; n = info{2}; c = info{3};
   for i=2:length(opList)
      info = opList{i}([],0);
      if (info{1} ~= m) || (info{2} ~= n)
         error(['Operator ' int2str(i) ' is not consistent with the first operators']);
      end

      c = c | info{3};
   end
   if any(~isreal(w)), c = ones(1,4); end;
end

op = @(x,mode) opSum_intrnl(m,n,c,opList,w,x,mode);


function y = opSum_intrnl(m,n,c,opList,w,x,mode)
checkDimensions(m,n,x,mode);
if mode == 0
   y = {m,n,c,{'Sum',opList}};
elseif mode == 1
   y = w(1) * opList{1}(x,mode);
   for i=2:length(opList)
     y = y + w(i)*opList{i}(x,mode);
   end
else
   y = conj(w(1)) * opList{1}(x,mode);
   for i=2:length(opList)
     y = y + conj(w(i))*opList{i}(x,mode);
   end
end
