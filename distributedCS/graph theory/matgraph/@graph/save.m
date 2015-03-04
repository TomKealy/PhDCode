function save(g,filename)
% save(g,filename) --- save a graph to disk
% The graph g is saved to a file named in the argument filename. 


fid = fopen(filename,'w');

if (fid == -1)
    error(['Cannot open "', filename, '" for output']);
end

n = nv(g);
m = ne(g);

fprintf(fid,'%%saved graph data\n');
fprintf(fid,'sp = %d;\n', issparse(g));
fprintf(fid,'nverts = %g;\n',n);
fprintf(fid,'nedges = %g;\n',m);

fprintf(fid,'a = zeros(%d,1);\n', m);
fprintf(fid,'b = zeros(%d,1);\n', m);

elist = edges(g);
a = elist(:,1);
b = elist(:,2);

step = 20;
for k=1:step:m
    last = min([k+step,m]);
    write_range(fid,a,k,last,'a');
    write_range(fid,b,k,last,'b');
end

fprintf(fid,'elist = [a,b];\n');

if hasxy(g)
    xy = getxy(g);
    x = xy(:,1);
    y = xy(:,2);
    fprintf(fid,'x = zeros(%d,1);\n', n);
    fprintf(fid,'y = zeros(%d,1);\n', n);
    
    step = 10;
    
    for k=1:step:n
        last = min([k+step,n]);
        write_range(fid,x,k,last,'x');
        write_range(fid,y,k,last,'y');
    end
    
    fprintf(fid,'xy = [x,y];\n');
else
    fprintf(fid,'xy = [];\n');
end

if is_labeled(g)
    fprintf(fid,'labs = {');
    for v=1:n
        fprintf(fid,['''', get_label(g,v), '''']);
        if (v<n)
            fprintf(fid,';');
        else
            fprintf(fid,'};\n');
        end
    end
else
    fprintf(fid,'labs={};\n');
end
fclose(fid);



function write_range(fid, array, first, last, name)

fprintf(fid,'%s(%d:%d) = [', name, first, last);
fprintf(fid,'%g;',array(first:(last-1)));
fprintf(fid,'%g];\n', array(last));
