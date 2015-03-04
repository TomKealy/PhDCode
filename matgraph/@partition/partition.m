function p = partition(param)
% partition --- constructor for the partition class
% p = partition(n) creates a default partition of [n] with all parts of
% size 1 (finest possible partition).
% p = partition(cell_array) creates a partition whose parts are given in
% the cell array.


if nargin==0
    param = 0;
end

% if called with a single argument, presumably a positive integer

if ~isa(param,'cell')
    n = param;
    p.array = logical(speye(n,n));  % singleton parts
    p.array = sortrows(p.array);
    p = class(p,'partition');
    return
end

% otherwise, called with a cell array containing the parts of p

% special case: param is an empty cell array
if isempty(param)
    p.array = speye(0);
    p = class(p,'partition');
    return
end

% Find nv (n) and np (m)
maxv = 0;
for i=1:length(param)
    m = max(param{i});
    maxv = max([m,maxv]);
end
n = maxv;
m = length(param);


% allocate the matrix
p.array = logical(sparse([],[],[],m,n,n));

% load the entries
for i=1:m
    row = zeros(1,n);
    row(param{i}) = 1;
    p.array(i,:) = logical(row);
end
p.array = sortrows(p.array);

p = class(p,'partition');

if ~check(p)
    p.array = speye(0);
    error('The cell array does not define a valid partition');
end

