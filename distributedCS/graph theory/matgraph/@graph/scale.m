function scale(g,s)
% scale(g,s) --- rescale the embedding of g by s

if ~hasxy(g)
    embed(g)
end

xy = getxy(g);
xy = s*xy;
embed(g,xy);