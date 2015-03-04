function ldraw(g,line_style)
% ldraw(g,line_style) --- draw a graph with vertices marked with their labels
% If the graph is unlabled, we use the vertex numbers instead.
% See also draw, cdraw, and ndraw.


if ~is_labeled(g)
    if nargin ==1
        ndraw(g);
        return
    else
        ndraw(g,line_style)
        return
    end
end

if nargin == 1
    draw(g);
else
    draw(g,line_style);
end

draw_labels(g);
