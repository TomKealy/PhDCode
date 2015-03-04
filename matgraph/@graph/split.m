function part = split(g,how)
% part = split(g) --- partition the vertices of g into two subsets
% The graph is partitioned based on the eigenvector associated with the 2nd
% smallest eigenvalue of the Laplacian of g. 
% split(g,'even') splits the graph into two equal (or as close as possible)
% parts

L = laplacian(g);

[V,d] = eig(L);
v2 = V(:,2);

a = 0; % place to cut v2

if (nargin == 2)
	if (strcmp(how, 'even'))
		a = mymedian(v2);
	end
end




part = cell(1,2);

part{1} = find(v2 > a);
part{2} = find(v2 <= a);
part = partition(part);

function m = mymedian(list) 

n = length(list);
if (mod(n,2)==1)
    m = median(list);
else
    list = sort(list);
    nn = n/2;
    m = (list(nn)+list(nn+1))/2;
end