function grid(g,a,b)
% grid(g,a,b) --- create an a-by-b grid graph
% grid(g,a) is the same as grid(g,a,a)



if nargin==2
    b = a;
end

n = a*b;

resize(g,0);
resize(g,n);

for i=1:a      % row index
    for j=1:b  % col index
        v = ij2v(i,j,a,b);
        
        % left 
        if (i>1)
            w = ij2v(i-1,j,a,b);
            add(g,v,w);
        end
        
        % right
        if (i<a)
            w = ij2v(i+1,j,a,b);
            add(g,v,w);
        end
        
        % down
        if (j>1)
            w = ij2v(i,j-1,a,b);
            add(g,v,w);
        end
        
        % up
        if (j<b)
            w = ij2v(i,j+1,a,b);
            add(g,v,w);
        end
    end
end


% now create the embedding

xy = zeros(n,2);

for i=1:a
    for j=1:b
        v = ij2v(i,j,a,b);
        xy(v,:) = [j,i];
    end
end
embed(g,xy)




% this is a helper function that converts i,j grid coordinates to a vertex
% number (1,1) --> 1 up to (a,b) --> ab

function v = ij2v(i,j,a,b)

v = i + a*(j-1);
