function icosahedron(g)
% icosahedron(g) --- overwrite g with the icosahedron graph

resize(g,12);
clear_edges(g);
full(g);

elist = [
    1   2
    1   3
    1   4
    1   5
    1   9
    2   3
    2   5
    2   6
    2   7
    3   7
    3   8
    3   9
    4   5
    4   9
    4   10
    4   12
    5   6
    5   10
    6   7
    6   10
    6   11
    7   8
    7   11
    8   9
    8   11
    8   12
    9   12
    10  11
    10  12
    11  12
    ];
add(g,elist);

t0 = 0;
t1 = -2*pi/3;
t2 = 2*t1;

outer = [
    sin(t0)     cos(t0)
    sin(t1)     cos(t1)
    sin(t2)     cos(t2)
    ];

hex = [];

for k=0:5
    hex = [hex; sin(-k*pi/3), cos(-k*pi/3)];
end

inner = outer*[cos(pi/3), sin(pi/3); -sin(pi/3), cos(pi/3)];


xy = [3*outer; hex; inner/3];
embed(g,xy)

    

