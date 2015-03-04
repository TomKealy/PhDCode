function k = plus(g,h)
% k = plus(g,h) --- disjoint union (sum) of graphs -- DEPROCATED!!
% Note: this creates a new graph k.
% May also be called simply as h = g+k.
% DANGER: Using this may result in a dangling handle and a memory leak. See
% the documentation. 
% THIS FUNCTION IS STRONGLY DEPROCATED. I PLAN TO REMOVE IT FROM MATGRAPH.

global GRAPH_MAGIC

k = graph;
disjoint_union(k,g,h);