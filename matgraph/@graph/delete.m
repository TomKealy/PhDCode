function delete(g,i,j)
% delete --- delete vertices or edges from a graph
% delete(g,v) --- delete vertex v
% delete(g,vlist) --- delete vertices in column vector v
% delete(g,i,j) --- delete edge ij
% delete(g,elist) --- delete edges in m-by-2 array elist

global GRAPH_MAGIC

n = nv(g);

if nargin==3   % simple delete an edge
   if (i<1) || (i>n) || (j<1) || (j>n) 
       return
   end
   GRAPH_MAGIC.graphs{g.idx}.array(i,j) = 0;
   GRAPH_MAGIC.graphs{g.idx}.array(j,i) = 0;
   return
end

% from here on, we're in the single argument case

[nr,nc] = size(i);

% delete a set of vertices (or a single vertex)

if  (nc==1)
    kill_list = i(find((i>0) & (i<=n)));
    keep_list = setdiff(1:n,kill_list);
    if hasxy(g)
        GRAPH_MAGIC.graphs{g.idx}.xy = ...
            GRAPH_MAGIC.graphs{g.idx}.xy(keep_list,:);
    end
    if is_labeled(g)
        GRAPH_MAGIC.graphs{g.idx}.labels = ...
            GRAPH_MAGIC.graphs{g.idx}.labels(keep_list);
    end
    GRAPH_MAGIC.graphs{g.idx}.array = ...
        GRAPH_MAGIC.graphs{g.idx}.array(keep_list,keep_list);
    return
end


if (nr==1) && (nc ~=2)
    error('for delete(g,vlist), vlist must be a single COLUMN vector');
end

% finally, we're in the delete a list of edges case

if (nc ~= 2)
    error('delete(elist) requires elist to have exactly two columns')
end

for k=1:nr
    delete(g,i(k,1),i(k,2));
end
    
   
