function resize(g, n)
% resize(g,n) --- change the number of vertices in g to n
% note: if n is less than nv(g), vertices will be lost
% if the graph has an embedding, all new vertices are sited at the origin. 


old_n = nv(g);

if (n == old_n)   % no change 
    return
end

global GRAPH_MAGIC

if (n < old_n)    % shrink
    GRAPH_MAGIC.graphs{g.idx}.array = ...
        GRAPH_MAGIC.graphs{g.idx}.array(1:n,1:n);
    make_logical(g);
    if hasxy(g)
        GRAPH_MAGIC.graphs{g.idx}.xy = ...
            GRAPH_MAGIC.graphs{g.idx}.xy(1:n,:);
    end
    if is_labeled(g)
        GRAPH_MAGIC.graphs{g.idx}.labels = ...
            GRAPH_MAGIC.graphs{g.idx}.labels(1:n);
    end
    return
end

% we need to add extra zeros to g's adjacency matrix. we decide what to do
% based on g's storage class.

if (isfull(g))
    A = zeros(old_n, n-old_n);
    B = zeros(n-old_n, n-old_n);
else
    A = sparse([],[],[],old_n,n-old_n);
    B = sparse([],[],[],n-old_n,n-old_n);
end

if hasxy(g)
    morexy = zeros(n-old_n,2);
    GRAPH_MAGIC.graphs{g.idx}.xy = ...
        [GRAPH_MAGIC.graphs{g.idx}.xy; morexy];
end
 
if is_labeled(g)
    newlabels = cell(n,1);
    for k=1:old_n
        newlabels{k} = GRAPH_MAGIC.graphs{g.idx}.labels{k};
    end
    for k=(old_n+1):n
        newlabels{k} = int2str(k);
    end
    GRAPH_MAGIC.graphs{g.idx}.labels = newlabels;
end

GRAPH_MAGIC.graphs{g.idx}.array = ...    
    [GRAPH_MAGIC.graphs{g.idx}.array, A; A', B];
make_logical(g);