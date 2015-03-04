function draw(g,line_style)
% draw(g) --- draw g in a figure window
% draw(g,line_style) --- lines have given line_style
% see also ndraw, ldraw, and cdraw
  
% edit these to change the colors 
edge_color = 'b';
vertex_color = 'r';
vertex_fill = 'w';
r = 0.15;
  
n = nv(g);
  
if nargin < 2
    line_style='-';
end


if ~hasxy(g)
    embed(g);
end

xy = getxy(g);

% first draw the edges

elist = edges(g);
for j=1:ne(g)
    u = elist(j,1);
    v = elist(j,2);
    x = xy([u,v],1);
    y = xy([u,v],2);
    line(x,y,'Color', edge_color,'LineStyle',line_style);
end


% now draw the vertices
  
for v=1:n
    x = xy(v,1);
    y = xy(v,2);
    rectangle('Position', [x-r/2, y-r/2, r, r],...
              'Curvature', [1 1], ...
              'EdgeColor', vertex_color, ...
              'FaceColor', vertex_fill);    
end



axis equal
axis off
