function dodecahedron(g)
% dodecahedron(g) --- overwrite g with the dodecahedron graph

resize(g,20)
clear_edges(g)
full(g)

elist = [
    1   2
    1   5
    1   6
    2   3
    2   7
    3   4
    3   8
    4   5
    4   9
    5   10
    6   11
    6   15
    7   11
    7   12
    8   12
    8   13
    9   13
    9   14
    10  14
    10  15
    11  16
    12  17
    13  18
    14  19
    15  20
    16  17
    16  20
    17  18
    18  19
    19  20
    ];

add(g,elist);

t0 = 0;
t1 = -2*pi/5;
t2 = 2*t1;
t3 = 3*t1;
t4 = 4*t1;

ring1 = [
    sin(t0)     cos(t0)
    sin(t1)     cos(t1)
    sin(t2)     cos(t2)
    sin(t3)     cos(t3)
    sin(t4)     cos(t4)
    ];

ring2 = ring1 * [cos(t1/2),-sin(t1/2);sin(t1/2),cos(t1/2)];

xy = [
    4*ring1;
    3*ring1;
    2*ring2;
    ring2
    ];

embed(g,xy)
