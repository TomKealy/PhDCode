function c = edge_color(g,algo)
% edge_color(g,algo) --- find a proper edge coloring of the graph g
% The returned matrix has three columns. Each row of this matrix of the
% form [i j c] where ij is an edge of the graph and c is a color number
% (from 1 to the number of colors used). 
% 
% algo denotes an algorithm. The choices are:
%
% greedy:    use a greedy coloring, based on the line graph (default)
% 
% optimal:   create an optimal coloring, based on the line graph
%
% vizing:    create a coloring that uses at most max(deg(g))+1 colors
%            NOTE: THIS IS UNIMPLEMENTED AT THIS TIME!

if nargin < 2
	algo = 'greedy';
end

switch (lower(algo))
	case {'greedy','optimal'}
		h = graph;
		line_graph(h,g);
		p = color(h,algo);
		c = [edges(g), array(p)'];
		free(h);
		
	case 'vizing'
		warning('edge_color:novizing',...
			'The "vizing" algorithm is not implemented. Using greedy algorithm.')
		c = edge_color(g,'greedy');
		
	otherwise
		error(['Unknown algorithm: ', algo]);
end