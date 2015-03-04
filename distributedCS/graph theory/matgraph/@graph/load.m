function load(g,filename)
% load(g,filename) --- read a saved graph on disk
% reads in a graph saved by a previous call to save(g,filename)


fid = fopen(filename,'r');

if fid == -1
    error(['File "', filename, '" cannot be opened for reading']);
end

A = textscan(fid,'%s\n','whitespace','\n');
A = A{1};

nl = length(A);

for k=1:nl
    line = A{k};
    eval(line);
end

resize(g,0);

if (sp == 1)
    sparse(g)
else
    full(g)
end

resize(g,nverts);
add(g,elist);

if (length(xy) > 0)
    embed(g,xy);
end

n = nv(g);

if (length(labs)>0)
    for k=1:n
        label(g,k,labs{k});
    end
end

fclose(fid);
