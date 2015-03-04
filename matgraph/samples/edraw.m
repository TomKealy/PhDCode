function edraw(g,elist)
n = nv(g);
h = graph(n);       % create a graph with same number of vertices as g
add(h,elist);       % add edges in elist to h
embed(h, getxy(g)); % copy g's embedding
draw(g,':');        % draw g with dotted lines
draw(h);            % overdraw h with solid lines
free(h);            % release h