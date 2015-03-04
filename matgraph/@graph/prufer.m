function output = prufer(g, code)
% prufer --- convert a tree to/from its Prufer code
% output = prufer(g) returns the Prufer code for g (assuming g is a tree)
% prufer(g,code) overwrites g with a tree based on the code 
%
% The Prufer code is a way to map bijectively trees on n vertices 
% into n-2 long sequences of integers drawn from [n].

if nargin==1
    n = nv(g);
    m = ne(g);
    if (m ~= n-1)
        error('Input graph is not a tree')
    end
    if (~isconnected(g))
        error('Input graph is not a tree')
    end
    
    if (n < 2)
        error('Algorithm applies only to trees with 2 or more vertices')
    end
    
    if (n==2)
        output=[];
        return
    end
    
    t = graph;  % temporary copy of g
    copy(t,g);
    
    names = 1:n;  % names of the vertices
    output = zeros(1,n-2);
    
    for i=1:n-2
        % find least leaf
        d = deg(t);
        v = min(find(d==1));
        
        % find neighbor of the least leaf
        w = neighbors(t,v);
        
        output(i) = names(w); % add this vertex to the output
        
        % delete v and its name
        delete(t,v);
        names = [names(1:v-1),names(v+1:end)];
    end
        
    free(t)
    return
end

if ~prufer_check(code)
    error('Input is not a valid Prufer code')
end

n = length(code) + 2;

resize(g,0);
resize(g,n);

verts = 1:n;
edges = zeros(n-1,2);

for k=1:n-2
    u = min(setdiff(verts,code));
    v = code(1);
    e = [u,v];
    edges(k,:) = e;
    code = code(2:end);
    i = find(verts == u);
    verts = [verts(1:i-1),verts(i+1:end)];
end


edges(n-1,:) = verts;
add(g,edges);
    



function ok = prufer_check(list)  % check if list is a valid Prufer code
ok = false;
n = length(list)+2;

if (any(list > n)) || (any(list < 1)) || (any(list ~= round(list)))
    return
end

ok = true;


