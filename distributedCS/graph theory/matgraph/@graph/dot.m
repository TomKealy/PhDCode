function dot(g,filename)
% dot(g,filename) --- save graph for use by graphviz.
% filename is the name of the file to save and should end with ".dot"
% writes the graph to disk in a form that can be read into the graphviz
% tools such as dot, neato, etc. 
% If the graph is labeled, we use the vertex labels instead of the vertex
% numbers. In this case, it is vital that the labels be distinct.

graph_name = ['Graph_', int2str(g.idx)];

if nargin==1
    filename = [graph_name,'.dot'];
end

fid = fopen(filename,'w');
if (fid<0)
    error(['Cannot open ', filename, ' for writing']);
end

label_state = is_labeled(g);
if (~label_state)
    label(g)
end

fprintf(fid,['strict graph ', graph_name, '{ \n']);

n = nv(g);

for v=1:n
    fprintf(fid,'"%s";\n', get_label(g,v));
end

elist = edges(g);
m = ne(g);

for k=1:m
    a = get_label(g,elist(k,1));
    b = get_label(g,elist(k,2));
    fprintf(fid, '"%s" -- "%s"; \n', a, b);
end

fprintf(fid,'}\n');


if (~label_state)
    clear_labels(g);
end

disp(['Wrote "', filename,'"']);
fclose(fid);