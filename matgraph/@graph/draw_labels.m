function draw_labels(g,offset)
% draw_labels(g) --- add labels to a drawing of g
% assuming g has already been drawn in a figure window (with draw or cdraw)
% add labels to the drawing. 
% Note: draw(g); add_labels(g) is tantamount to ldraw(g)
% If the drawing already has vertex numbers, this might make a mess.
%
% An optional second argument gives an offset of the label away from the
% vertex: draw_labels(g, [dx, dy]) moves the label dx to the right and dy
% up. Suggested offset: [0.1, -0.1]

if (~is_labeled(g))
    label(g)
end

if nargin == 1
    offset = [0 0];
end

xy = getxy(g);
n = nv(g);

for v=1:n
    x = xy(v,1) + offset(1);
    y = xy(v,2) + offset(2);
    text(x,y,get_label(g,v)); 
end