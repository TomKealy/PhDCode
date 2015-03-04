function yn = isfull(g)
% isfull(g) --- check if g's adjacency matrix is full
yn = ~issparse(g);