function blist = bridges(g,algo)
% bridges(g,algo) --- find all cut edges in g
% algo specifies the algorithm. Current choices are these:
%
% 'path'    Deletes each edge from the graph and checks if there is a path
%           between its endpoints. A few tricks are employed so we don't
%           have to check every edge.
%
% 'matrix'  Uses a basis for the null space of the signed incidence
%           matrix. 
%
% No decision yet on which of these is better as the default. 
%
% Main idea for the path code due to Mel Janowitz.

DEFAULT_ALGORITHM = 'matrix';

if nargin < 2
	algo = DEFAULT_ALGORITHM;
end

switch lower(algo)
	case {'path', 'paths'}
		blist = path_bridges(g);
	case 'matrix'
		blist = matrix_bridges(g);
	otherwise
		error(['Unknown algorithm: ', algo]);
end

function blist = matrix_bridges(g)
% Based on the following observation. Let M be the signed incidence matrix
% of the graph. An edge e is a bridge of g if and only if every vector in
% the null space of M has a 0 in the e'th coordinate. 

M = full(incidence_matrix(g,'signed'));
N = null(M);

% find rows that are entirely zero (or nearly so)
N = abs(N) < 100*eps;
idx = all(N,2);

blist = edges(g);
blist = blist(idx,:);


function blist = path_bridges(g)
% Implement the path method for finding bridges

bgraph = graph;     % hold the cut edges
notbridges = graph; % hold the non cut edges
gsave = graph;      % save a copy of g since our code modifies g
copy(gsave,g);

elist = edges(g);  % list of edges
m = ne(g);         % number of edges

for k=1:m
	% xy is the edge we're considering
	x = elist(k,1);
	y = elist(k,2);
	
	% if xy is in a cycle, we skip the rest of this iteration
	if (has(notbridges,x,y))
		continue
	end
	
	% find an xy-path in G-e
	delete(g,x,y);
	P = find_path(g,x,y);
	
	
	% if no such path, we have a cut edge
	if isempty(P)
		add(bgraph,x,y)
	else
		% otherwise, note that xy and all edges on P cannot be cut edges
		add(g,x,y);
		add(notbridges,x,y)
		for j=1:length(P)-1
			add(notbridges,P(j),P(j+1))
		end
	end
end

blist = edges(bgraph);

% restore g
copy(g, gsave)
% release temporary storage
free(bgraph)
free(notbridges)
free(gsave)
