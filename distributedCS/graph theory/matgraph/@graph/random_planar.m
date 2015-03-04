function random_planar(g,n)
% random_planar(g,n) --- create a random planar triangulation
% This does NOT create a random "uniform" triangulation. Rather, it places
% n points at random in the plane and builds the Delaunay triangulation
% thereon. 

xy = rand(n,2)*sqrt(n);
resize(g,0)
resize(g,n)
embed(g,xy)

triangles = delaunay(xy(:,1), xy(:,2));
nt = size(triangles,1);

edges = [triangles(:,1:2); triangles(:,2:3); triangles(:,[1,3])];
edges = sortrows(edges);
edges = unique(edges,'rows');
add(g,edges)