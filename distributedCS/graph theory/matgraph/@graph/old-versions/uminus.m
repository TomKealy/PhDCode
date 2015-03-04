function h = uminus(g)
% uminus(g) --- operator form of complement: h = -g;  DEPROCATED!!!!
% WARNING: creates a new graph
% THIS FUNCTION IS STRONGLY DEPROCATED. I PLAN TO REMOVE IT FROM MATGRAPH.

h = graph;
copy(h,g);
complement(h);