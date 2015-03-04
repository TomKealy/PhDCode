function add(g,i,j)
% add --- add edge(s) to the graph
% add(g,i,j) --- adds the edge ij
% add(g,elist) --- adds edges specified in elist to the graph
%                  elist is an n-by-2 array of edges
% If an end point is larger than the number of vertices, the graph is
% extended. End points must be positive integers. Loops are ignored.

global GRAPH_MAGIC
n = nv(g);


if nargin == 3
    if (i==j) || (i<1) || (j<1) 
        return
    end
    
    maxij = max(i,j);
    if maxij > n
        resize(g,maxij)
    end

    GRAPH_MAGIC.graphs{g.idx}.array(i,j) = 1;
    GRAPH_MAGIC.graphs{g.idx}.array(j,i) = 1;
    return
end

if nargin == 2
    [nr,nc] = size(i);
    % check shape
    if (nc ~= 2)
        error('in add(elist), elist should have exactly two columns');
    end
    
    maxv = max(max(i));
    if maxv > n
        if maxv > set_large
            sparse(g)
        end
        resize(g,maxv)
    end
    n = nv(g);
    was_full = isfull(g);
    
    data = [i; n,n]; % note phony last row
    m = size(data,1);
    data = [data, ones(m,1)];
    A = spconvert(double(data));
    for k=1:n
        A(k,k)=0;
    end
    
    A = (A+A')>0;
    A = logical(A);
    if (was_full) 
        A = full(A);
    end
    GRAPH_MAGIC.graphs{g.idx}.array = GRAPH_MAGIC.graphs{g.idx}.array | A;
    return
end
