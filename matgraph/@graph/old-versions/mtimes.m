function h = mtimes(g1,g2)
% h = mtimes(g1,g2) [or g1*g2] --- operator form of cartesian DEPROCATED!!
% WARNING: This method creates a new graph; if you don't save the result, a
% dangling pointer / memory leak will result.
% THIS FUNCTION IS STRONGLY DEPROCATED; I PLAN TO REMOVE IT FROM MATGRAPH.

h = graph;
cartesian(h,g1,g2);