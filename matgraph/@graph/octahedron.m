function octahedron(g)
% octahedron(g) --- overwrite g with the octahedron graph, K(2,2,2)


resize(g,6);
full(g);
clear_edges(g);

elist = [
    1   2
    1   3
    2   3
    4   5
    4   6
    5   6
    1   5  
    3   5
    1   6
    2   6
    2   4
    3   4
];

add(g,elist);


xy = [
    0               1/2
    -sqrt(3)/4       -1/4
    sqrt(3)/4      -1/4
    0               -8/3
    4/sqrt(3)       4/3
    -4/sqrt(3)      4/3
    ];

embed(g,xy);