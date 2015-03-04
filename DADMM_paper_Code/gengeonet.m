function [Network] = gengeonet(P, d)

fid = fopen('net_data.dat','r');
data = textscan(fid, '%s');
fclose(fid);
data = regexprep(data{1},'[\[\],]','');

ind = find(strcmp(data, 'conn'));
conn = str2double(reshape(data(ind+2),1,1)');

beg = find(strcmp(data, 'Adj'));
posind = find(strcmp(data, 'pos'));

Adj = str2double(reshape(data(beg+2:posind-1),50,50)');
pos = str2double(reshape(data(posind+2:end),2,50)');

parameters = cell(1,1);
parameters{1} = d;

Network = struct('P', P, 'Adj', Adj, 'pos', pos, 'Type', 'Geometric', 'Parameters', parameters);

P = Network.P;
Adj = Network.Adj;
Type = Network.Type;
Parameters = Network.Parameters;

E = [];

e = 1;

NAdj = Adj + Adj';

for i = 1 : P
    for j = i+1 : P
        if NAdj(i,j) > 0
            v = [i;j];
            E = [E , v]; %#ok<AGROW>
            e = e + 1;
        end
    end
end

% Determine the color of the graph
clear global GRAPH_MAGIC;
graph_init
g_matgraph = graph(P);
add(g_matgraph,E');
if sum((Adj - matrix(g_matgraph)) ~= zeros(size(Adj)))
    error('Adjacency matrices given by R and matgraph differ!');
end
partition_by_colors_aux = color(g_matgraph, 'greedy');
Diameter = diam(g_matgraph);
free(g_matgraph);
clear global GRAPH_MAGIC;

partition_by_colors = parts(partition_by_colors_aux);
Num_Colors = length(partition_by_colors);

Neighbors = cell(P,1);
Degrees = zeros(P,1);

for p = 1 : P
    Neighbors{p} = find(Adj(p,:));
    Degrees(p) = length(Neighbors{p});
end

Max_Degree = max(Degrees);

Network = struct('P', {P}, ...
    'Adj', {Adj}, ...
    'Type', {Type}, ...
    'Parameters', {Parameters}, ...
    'Num_Colors', {Num_Colors}, ...
    'Partition', {partition_by_colors}, ...
    'Diameter', {Diameter}, ...
    'Neighbors', {Neighbors}, ...
    'Degrees', {Degrees}, ...
    'Max_Degree', {Max_Degree}, ...
    'Positions', {pos} ...
    );
end