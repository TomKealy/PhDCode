function make_logical(g)
% make_logical(g) --- ensure that the internal storage for g's matrix is a
% logical matrix (for space efficiency). 
% This function is called by various matgraph routines and it should not 
% be necessary to call this in ordinary use of matgraph.

global GRAPH_MAGIC

if (~islogical(GRAPH_MAGIC.graphs{g.idx}.array))
    GRAPH_MAGIC.graphs{g.idx}.array = ...
        logical(GRAPH_MAGIC.graphs{g.idx}.array);
end