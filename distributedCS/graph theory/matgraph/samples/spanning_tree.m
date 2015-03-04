%% Finding spanning trees
% Matgraph provides two functions for finding spanning trees:
% 
% * |bfsptree| --- breadth first spanning tree
% * |dfsptree| --- depth first spanning tree
%% Set up the graph
% We illustrate finding spanning trees in the icosahedron graph.
g = graph;
icosahedron(g);
%% Depth first spanning tree
t = graph;
dfstree(t,g);
clf; draw(g,':');
draw(t)
%% Breadth first spanning tree
bfstree(t,g);
clf; draw(g,':');
draw(t)
%% Counting spanning trees
% The matrix-tree theorem makes enumerating spanning trees in Matlab
% especially easy.
nt = nsptrees(g);
disp(['The icosahedron graph has ', int2str(nt),' spanning trees']);
%% Release the graphs
free(t)
free(g)