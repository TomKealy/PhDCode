function idx = find_row(r, M, a, b)
% find_row(r,M)
% find row r in matrix M. Assume rows of M are sorted and unique
%
% find_row(r,M,a,b)
% a = lowest row to check (default 1)
% b = highest row to check (default number of rows in M)
%
% return the idx of the row, or 0 if not found

if nargin < 3
	a = 1;
	b = size(M,1);
end

if a>b
	idx = 0;
	return
end

if a==b
	if vector_compare(r,M(a,:)) == 0
		idx = a;
	else
		idx = 0;
	end
	return
end

mid = floor( (a+b)/2 );

sg = vector_compare(r, M(mid,:));

switch sg
	case 0,
		idx = mid;
		return;
	case -1,
		idx = find_row(r, M, a, mid-1);
		return
	case 1,
		idx = find_row(r, M, mid+1, b);
end
