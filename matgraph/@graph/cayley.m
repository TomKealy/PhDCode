function cayley(g,perms,verbose)
% cayley(g,perms) -- create a Cayley graph (undirected)
% g is the graph to be written
% perms is a cell array of permutations that are the generators of a group.
% The vertices of g are the elements of the generated group. There is an
% edge from u to v in g provided there is a generator x so that u=v*x or 
% v = u*x.
%
% cayley(g,perms,verbose) with verbose = true displays progress as this
% function works. [Default is verbose = false.]

if nargin<3
	verbose = false;
end

if ~isa(perms,'cell')
	error('2nd argument must be a cell array of permutations');
end

perms = perms(:);    % make it 1-dimensional
np = length(perms);  % number of perms

% if there are no perms in the list, then the group is {1} and the graph is
% a single vertex
if np == 0
	resize(g,1);
	return
end

if verbose
	disp('Generating the group elements from these:')
end

group = gen_group(perms,verbose);

if verbose
	disp('Generating edge list')
end

edges = [];

ng = size(group,1);  % size of the group

for i=1:ng
	u = group(i,:);
	pu = permutation(u);
	for k=1:np
		pv = pu * perms{k};
		v = array(pv);
		j = find_row(v,group);
		edges = [edges; [i,j]];
	end
end

edges = unique(sortrows(edges),'rows');

m = size(edges,1);
if verbose
	disp([int2str(m),' edges created']);
end

resize(g,0)
add(g,edges);

% label the nodes
clear_labels(g);
for u = 1:nv(g)
	str = display(permutation(group(u,:)));
	label(g,u,str);
end


%%
function group = gen_group(perms,verbose)
% create the full list of elements of a group based on a matrix of
% generators.


np = length(perms);
psize = 0;

for i=1:np
	pis = length(perms{i});
	if pis > psize
		psize = pis;
	end
end


gen_matrix = zeros(np, psize);
for i=1:np
	p = perms{i} * permutation(psize);
	perms{i} = p;
	row = array(p);
	gen_matrix(i,:) = row;
	if verbose
		display(perms{i})
	end
end

% initially, the group is the identity and all the generators
group = [1:psize; gen_matrix];
group = sortrows(group);
last = group;

while ~isempty(last)
	ng = size(group,1);
	if verbose
		disp(['Expanding. Group size = ', int2str(ng)]);
	end
	new_rows = [];
	for k=1:size(last,1)
		p = permutation(last(k,:));
		for j=1:np
			q = p * perms{j};
			row = array(q);
			idx = find_row(row,group);
			if idx == 0
				new_rows = [new_rows;row];
			end
		end
	end
	if size(new_rows,1) == 0
		return
	end
		
	new_rows = unique(sortrows(new_rows),'rows');
	last = [];
	for k=1:size(new_rows,1)
		r = new_rows(k,:);
		if find_row(r,group) == 0
			last = [last;r];
		end
	end
	group = [group; last];
	group = unique(sortrows(group),'rows');
	
end

