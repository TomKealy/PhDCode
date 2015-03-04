function mobius(g,nverts)
% mobius(g,nverts) --- create a Mobius ladder graph
% g is the graph to be written
% nverts is the number of vertices (must be even)


if (mod(nverts,2)==1)
	error('Number of vertices must be even')
end
if (nverts < 6)
	error('Number of vertices must be at least 6')
end

n = nverts/2;
cycle(g,nverts)
more_edges = zeros(n,2);
for v=1:n
	more_edges(v,:) = [v, v+n];
end

add(g,more_edges);

% create an embedding
rmxy(g);

% inner ring of vertices
t = (0:n-1)*2*pi/n;
x = n * cos(t)/6;
y = n * sin(t)/6;
xy1 = [x',y'];

% outer ring of vertices
x = (n+6)*cos(t)/6;
y = (n+6)*sin(t)/6;
xy2 = [x',y'];

xy = [xy1;xy2];

% rotate slightly so crossing is symmetric across the x-axis
theta = pi/n;
ct = cos(theta);
st = sin(theta);
rot = [ct,st; -st, ct];
xy = xy*rot;

embed(g,xy);
