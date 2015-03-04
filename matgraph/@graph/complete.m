function complete(g,a,b)
% complete: build complete and complete multipartite graphs
% complete(g) --- convert g to a complete graph on same vertex set
% complete(g,a) --- convert g to a complete graph on a vertices
% complete(g,a,b) --- convert g to K(a,b)
% complete(g,list) --- convert g to a multipartite graph

global GRAPH_MAGIC

% convert existing graph to complete

if nargin==1
    n = nv(g);
    fast_set_matrix(g,(ones(n) - eye(n)));
    return
end

% overwrite with K_a

if (nargin==2) && (length(a)==1)
    fast_set_matrix(g,ones(a)-eye(a));
    return
end

% overwrite with K(a,b)

if (nargin==3)
    A = zeros(a,a);
    B = zeros(b,b);
    Z = ones(a,b);
    fast_set_matrix(g,[A,Z;Z',B]);
    return
end

% last case: complete multipartite graph (a is a list)

n = sum(a);

GRAPH_MAGIC.graphs{g.idx}.array = true(n);
a = a(:)';
aa = [0,cumsum(a)];
for k=1:length(a)
    GRAPH_MAGIC.graphs{g.idx}.array(aa(k)+1:aa(k+1),aa(k)+1:aa(k+1))=0;
end


