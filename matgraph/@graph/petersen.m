function petersen(g)
% petersen(g) --- overwrite g with the Petersen graph

resize(g,0);
resize(g,10);
full(g);

elist = [
    1 6
    2 7
    3 8
    4 9
    5 10
    1 2
    2 3
    3 4
    4 5
    5 1
    6 8
    7 9
    8 10
    9 6
    10 7
    ];
add(g,elist);

t = (0:4)*2*pi/5;
x = sin(t);
y = cos(t);
xy1 = [x', y'];

xy = [2*xy1; xy1];
embed(g,xy);
