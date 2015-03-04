function random_regular(g,n,k)
% random_regular(g,n,k) --- create a random regular graph
% g is the graph to be overwritten
% n is the number of vertices
% k is the degree
%
% Note: This is a proof-of-concept version; no claims that this
% implementation is efficient! 

if mod(n*k,2)==1
    error('We number have n*k even or no regular graph is possible');
end

if (k>n-1)
    error('Degree too large');
end

if (k<0)
    error('Degree must be nonnegative')
end

h = graph(n);  % scratch graph

while(true)
    clear_edges(h);
    p = randperm(n*k);
    elist = reshape(p,n*k/2,2);
    elist = mod(elist,n)+1;
    if any(elist(:,1)==elist(:,2))  % reject if there's a loop
        continue
    end
    add(h,elist);
    if ne(h)== n*k/2
        break
    end
end

copy(g,h);
rmxy(g);
clear_labels(g);

free(h);