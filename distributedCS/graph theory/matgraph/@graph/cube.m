function cube(g,k)
% cube(g,k) --- create a k-cube (default k = 3)

if nargin==1
    cube(g,3);
    xy = [ 
        0   0
        1   0
        0   1  
        1   1
        ];
    xy = [xy ; 2*xy-1/2];
    embed(g,xy);
    return
        
end

if (k==1)
    path(g,2);
    return
end

if (k==2)
    resize(g,4);
    clear_edges(g);
    full(g);
    elist = [ 1 2; 2 4; 4 3; 3 1];
    add(g,elist);
    xy = [ 0 0 ; 1 0 ; 0 1 ; 1 1];
    embed(g,xy);
    return
end

cube(g,k-1);

k2 = graph(2);
tmp = graph;

complete(k2);
embed(k2,[-1,0;1,0]);
copy(tmp,g);
cartesian(g,tmp,k2);

free(k2);
free(tmp);
