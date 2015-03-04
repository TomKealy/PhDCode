function ndraw(g,line_style)
% ndraw(g) --- draw g in a figure window with numbered vertices
% ndraw(g,line_style) --- lines have given line_style
% see also draw, ldraw, and cdraw

if nargin == 1
    draw(g);
else
    draw(g,line_style);
end

xy = getxy(g);
n = nv(g);

for v=1:n
    x = xy(v,1);
    y = xy(v,2);
    text(x-.05,y,int2str(v)); 
end