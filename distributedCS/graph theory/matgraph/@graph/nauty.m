function nauty(g, filename)
% nauty(g,filename) -- save a graph in nauty format
% g is the graph, filename is a string containing the filename.
% It is recommended that the file name end ".dre"
%
% nauty is a program for computing automorphisms of graphs. It can be
% downloaded from http://cs.anu.edu.au/~bdm/nauty/
%
% Inside the dreadnaut program, type <filename to load the graph saved by
% this function.


fid = fopen(filename,'w');

if (fid == -1)
    error(['Cannot open "', filename, '" for output']);
end

n = nv(g);

fprintf(fid,'$ 1\n');     % sets nauty to start number vertices from 1
fprintf(fid,'n=%d\n', n); % write the number of vertices
fprintf(fid,'g\n');       % begin the graph

for v=1:n
	Nv = neighbors(g,v);
	Nv = Nv(Nv>v);
	if isempty(Nv)
		continue
	end
	fprintf(fid,'%d: ', v);
	fprintf(fid,'%d ', Nv);
	fprintf(fid,';\n');
end


fprintf(fid,'.\n');       % end writing the graph

