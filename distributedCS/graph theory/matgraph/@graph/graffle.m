function graffle(g, filename, width, rad)
% graffle(g, filename, width, rad) --- write graph in OmniGraffle format
%
% This writes a graph to the disk with the file name specified in the
% second argument (default is 'graph.graffle'). 
% The width (3rd argument) gives the overall size of the plot
% (default=450).
% The radius (4th argument) gives the size of a vertex (default=12).
% Note units are points = 1/72 inch. So 18 = 1/4 inch.
% 
% The output file can then be opened in OmniGraffle, a graph drawing 
% program for Macintosh. 

DEFAULT_FILE_NAME = 'graph.graffle';
DEFAULT_RADIUS = 12;
PAGE_SIZE = 450;

if nargin == 1
    filename = DEFAULT_FILE_NAME;
end

if nargin < 3
    width = PAGE_SIZE;
end

if nargin < 4
    rad = DEFAULT_RADIUS;
end

fid = fopen(filename,'w');

if (fid == -1)
    error(['Cannot open ', filename, ' for output']);
end

if (~hasxy(g))
    embed(g)
end

xy = getxy(g);
x = xy(:,1);
y = xy(:,2);


% move lower left to 0,0
x = x-min(x);
y = -y;
y = y-min(y);

% rescale to fit in reasonable square
m = max([x;y]);
scale = width/m;

x = round(scale*x)+5;
y = round(scale*y)+5;


% write XML header line
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');

% write start of XML GraphicsList array
fprintf(fid,'<dict>\n');
write_key(fid, 'GraphicsList');
fprintf(fid,'<array>');


% write vertices to XML file
for v = 1:nv(g)
    write_vertex(fid, v, x(v), y(v), rad)
end

% write edges
elist = edges(g);
m = ne(g);
for i=1:m
    u = elist(i,1);
    v = elist(i,2);
    write_edge(fid, u, v);
end


% end the array
fprintf(fid,'</array>\n');
% end the dictionary
fprintf(fid,'</dict>\n');

fclose(fid);



% page area: 
% x values run 0 to 520 or so
% y values run 0 to 700 or so


function write_key(fid,key_name)
fprintf(fid,[' <key>', key_name, '</key>\n']);

function write_string(fid, string)
fprintf(fid,[' <string>', string, '</string>\n']);

function write_integer(fid,n)
fprintf(fid,[' <integer>', int2str(n), '</integer>\n']);


% write vertex v @ coordinates (x,y) with radius r to output
function write_vertex(fid, v, x, y, r)
fprintf(fid,'<dict>\n');
write_key(fid,'Bounds');
fprintf(fid,[' <string>{{', int2str(x), ',', int2str(y), '}, {', ...
    int2str(r), ',', int2str(r),'}}</string>\n']);
write_key(fid,'Class');
write_string(fid, 'ShapedGraphic');
write_key(fid,'ID');
write_integer(fid,v);
write_key(fid,'Shape');
write_string(fid,'Circle');
write_key(fid,'Style');
fprintf(fid,' <dict>\n');
write_key(fid,'shadow');
fprintf(fid,' <dict>\n');
write_key(fid,'Draws');
write_string(fid, 'NO');
fprintf(fid,' </dict>\n </dict>\n</dict>\n\n');




% write an edge
function write_edge(fid, v, w)
fprintf(fid,'<dict>\n');
write_key(fid,'Class');
write_string(fid,'LineGraphic');

write_key(fid,'Head');
fprintf(fid,' <dict>\n');
write_key(fid,'ID');
write_integer(fid,v);
fprintf(fid,' </dict>\n');


write_key(fid,'Tail');
fprintf(fid,' <dict>\n');
write_key(fid,'ID');
write_integer(fid,w);
fprintf(fid,' </dict>\n');

fprintf(fid,'</dict>\n');


