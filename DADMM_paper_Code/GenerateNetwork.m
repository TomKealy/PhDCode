function [Network] = GenerateNetwork(num_nodes, type, parameters)

% [Network] = GenerateNetwork(num_nodes, type, parameters)
%
% Generates a connected undirected network with P = num_nodes. This 
% function works on:
%   - Windows: using an R interface (the graphs are generated with R)
%   - Linux: using a python package called Networkx (the graphs are 
%            generated in Networkx) 
%
% Since the input parameters for the network (explained ahead) might lead
% to unconnected networks, this function might change them in order to make
% the network connected with larger probability. Don't worry: you will know
% what the parameters that generated the network were because they will be
% stored in the output struct Network.
%
% This function requires:
%
%   matgraph: available at http://www.ams.jhu.edu/~ers/matgraph/
%
%   In Windows: R, R package igraph, and R-(D)Com interface (see README.txt
%               for installation)
%
%   In Linux: Python and NetworkX (see README.txt for installation)
%
%   
%
% Inputs:
% 
%   - num_nodes: number of nodes, or shortly, P
%   - type: (string) type of the network; possible types:
%           . Erdos-Renyi
%           . Watts-Strogatz
%           . Barabasi
%           . Geometric
%           . Lattice
%
%   - parameters: it is a cell with the required parameters to create the
%                 network; this depends on the type of network (see below)
%
%
% Outputs:
%  
%   - Network: is a struct with the following fields:
%           . P: number of nodes
%           . Adj: adjacency matrix: P x P
%           . Type: string with the name/type of the network
%           . Parameters: cell with the parameters of the network (note:
%                         the input parameters may be changed to create a
%                         connected network)
%           . Num_Colors: number of colors used to proper color the graph
%                         (using the matgraph package for Matlab)
%           . Partition: cell Num_colors x 1 with the coloring scheme
%           . Diameter: (number) diameter of the network
%           . Neighbors: cell P x 1, where the pth entry contains a vector
%                        of the size of the number of neighbors of node p;
%                        each vector contains the neighbors of node p
%           . Degrees: Vector P x 1, with the degree of node p at entry p
%           . Max_Degree: (number) maximum degree over all nodes
%
%
% Types of networks:
%
%   - Erdos-Renyi(p): Given P nodes, it goes to every pair of nodes and
%                     connects them with probability p. It creates a
%                     connected network with very high probability for
%                     p>log(P)/P
%
%   - Watts-Strogatz(nei,p): Creates a "ring" where each neighbor is
%                            connected with its closest 2*nei neighbors.
%                            Then, every link is rewired with probability
%                            p. Note: nei has to be greater than 1.
%
%   - Barabasi(m): It creates a scale free network. Each node is connected 
%                  to the previous network. m is the number of nodes to
%                  each new node connects.
%
%   - Geometric(d): It creates a square [0,1]x[0,1] and places P nodes 
%                   randomly. Then is connects every two nodes if their
%                   distance is less or equal to d. If the network fails to
%                   be connected, the parameter d is increased and the
%                   experiment repeated.
%
%   - Lattice: Creates a rectangular lattice with P nodes. To determine the
%              number of nodes in each direction, we do the following (to
%              try to make the lattice as square as possible): starting at
%              u = floor(sqrt(P)), we decrease u by one unit until it
%              reaches 1. If there is any number for which the P is
%              divisible for u, we stop. The other dimension will be P/u.
%              If it is a prime number, we will have a line. It returns in
%              parameters a vector with [n1,n2], the dimensions of the
%              grid.
%
%
% Note: All types of networks, except the Lattice, are generated in R 
%       (Windows) or in Networkx (Python)
%       
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Please send any questions, comments, or bug reports to joaomota@cmu.edu.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



P = num_nodes;
if P < 2
    error('Number of nodes has to be greater than 1');
end


if strcmp(type,'Erdos-Renyi')
    
    % Quantity we add to p everytime the network is not connected
    add_p = 0.01;   
    
    MAX_ITER = 10; % Maximum number of iterations (when the network is not 
                   % connected)
    
    if isempty(parameters)
        error('parameters field is empty');
    end
    p = parameters{1};
    if p < log(P)/P
        fprintf('Setting the Erdos-Renyi parameter to 1.01*log(P)/P = %f\n', 1.01*log(P)/P);
        p = 1.01*log(P)/P;
    end
    
    if ispc == 1        % If OS is Windows: use Matlab R-link software
        [status,msg] = openR;
        if status ~= 1
            disp(['Problem connecting to R: ' msg]);
        end
    end
    
    
    for iter = 1 : MAX_ITER

        [conn, Adj] = Erdos_Renyi_interface(P,p);
        
        if conn == 1 
            break;
        else
            p = p + add_p;
        end

    end
    
    if ispc == 1     % If OS is Windows
        closeR
    end
    
    if conn == 0
        error('Generated network (Erdos-Renyi) was not connected');
    end
    
    
    parameters = cell(1,1);
    parameters{1} = p;

    Network = struct('P', P, 'Adj', Adj, 'Type', 'Erdos-Renyi', 'Parameters', parameters);
    
    
elseif strcmp(type,'Watts-Strogatz')
    
    % Quantity we add to nei everytime the network is not connected
    add_nei = 1;   
    
    MAX_ITER = 10; % Maximum number of iterations (when the network is not 
                   % connected)
    
    if isempty(parameters)
        error('parameters field is empty');
    end
    nei = parameters{1}(1);
    p = parameters{1}(2);
    
    if nei <= 1 || p > 1 || p < 0
        error('Please read the definition of the Watts-Strogatz model: the parameters are not well set');
    end
    
    if ispc == 1     % If OS is Windows
        [status,msg] = openR;
        if status ~= 1
            disp(['Problem connecting to R: ' msg]);
        end
    end

    for iter = 1 : MAX_ITER

        [conn, Adj] = Watts_Strogatz_interface(P, nei, p);
        
        if conn == 1 
            break;
        else
            nei = nei + add_nei;
        end

    end
    
    if ispc == 1     % If OS is Windows
        closeR
    end

    if conn == 0
        error('Generated network (Watts-Strogatz) was not connected');
    end

    parameters = cell(1,1);
    parameters{1} = [nei, p];
    
    Network = struct('P', P, 'Adj', Adj, 'Type', 'Watts-Strogatz', 'Parameters', parameters);
    
    
elseif strcmp(type,'Barabasi')
    
    if isempty(parameters)
        error('parameters field is empty');
    end
    m = parameters{1};
    if m < 1
        fprintf('Setting the Barabasi parameter to m = 1\n');
        m = 1;
    end
    
    if ispc == 1     % If OS is Windows
        [status,msg] = openR;
        if status ~= 1
            disp(['Problem connecting to R: ' msg]);
        end
    end
    
    [conn, Adj] = Barabasi_interface(P, m);
    
    if ispc == 1     % If OS is Windows
        closeR
    end
    
    if conn == 0
        error('Generated network (Barabasi) was not connected');
    end

    parameters = cell(1,1);
    parameters{1} = m;
    Network = struct('P', P, 'Adj', Adj, 'Type', 'Barabasi', 'Parameters', parameters);
        
    
elseif strcmp(type,'Geometric')
    
    % Quantity we add to d everytime the network is not connected
    add_d = 0.02;   
    
    MAX_ITER = 50; % Maximum number of iterations (when the network is not 
                   % connected)
    
    if isempty(parameters)
        error('parameters field is empty');
    end
    d = parameters{1};
    
    if d <= 0
        error('Parameter d of geometric model is less than 0');
    end
    
    if ispc == 1     % If OS is Windows
        [status,msg] = openR;
        if status ~= 1
            disp(['Problem connecting to R: ' msg]);
        end
    end

    for iter = 1 : MAX_ITER

        [conn, Adj] = Geometric_interface(P, d);
        
        if conn == 1 
            break;
        else
            if mod(iter,5) == 0
                d = d + add_d;
            end
        end

    end
    
    if ispc == 1     % If OS is Windows
        closeR
    end
    
    if conn == 0
        error('Generated network (geometric) was not connected');
    end


    parameters = cell(1,1);
    parameters{1} = d;
    
    Network = struct('P', P, 'Adj', Adj, 'Type', 'Geometric', 'Parameters', parameters);
    
elseif strcmp(type,'Lattice')
    
    u = floor(sqrt(P));
    
    for i = u : -1 : 1
        if mod(P,i)==0
            break;
        end
    end
    
    n1 = i;
    n2 = P/i;
    
    if P ~= n1*n2
        error('Weird error: P is different than n1*n2 in Lattice type of network');
    end
    
    if ispc == 1     % If OS is Windows
        [status,msg] = openR;
        if status ~= 1
            disp(['Problem connecting to R: ' msg]);
        end
    end
    
    [conn, Adj] = Lattice_interface(n1, n2);
    
    if ispc == 1     % If OS is Windows
        closeR
    end
    
    if conn == 0
        error('Generated network (Lattice) was not connected');
    end

    parameters = cell(1,1);
    parameters{1} = [n1,n2];
    Network = struct('P', P, 'Adj', Adj, 'Type', 'Lattice', 'Parameters', parameters);
            
    
    
else
    error('Network type nonexistent');
end

% =========================================================================
% Determine the coloring scheme of the network, plus the other statistics

P = Network.P;
Adj = Network.Adj;
Type = Network.Type;
Parameters = Network.Parameters;

E = GN_convertAdj2ListOfNodes(Adj, P);

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
    'Max_Degree', {Max_Degree} ...
    );

% =========================================================================

end


function [conn, Adj] = read_Adjacency_and_connectivity(filename, P)
% [conn, Adj] = read_Adjacency_and connectivity(filename)
%
% This function is used by all functions below as a part of the R
% interface, namely when the OS is linux. All functions below produce a
% file with name 'filename' that contains:
%
% conn = val
% Adj =
% val val val ... val
% val val val ... val
% ...
% val val val ... val
%
% where val is either 0 or 1. The size of Adj matrix is P x P. This
% function returns both conn and Adj as output.

Adj = zeros(P,P);

fid = fopen(filename, 'r');
data = textscan(fid, '%s');
fclose(fid);

conn = str2double(data{1}(3));

index = 5;

for i = 1 : P
    for j = 1 : P
        Adj(i,j) = str2double(data{1}(index + P*(i-1) + j));
    end
end

end



function [conn, Adj] = Erdos_Renyi_interface(P,p)

% Python filename (for Linux)
Erdos_Renyi_filename = 'Erdos_Renyi.py';

% Filename of output file in case we execute the previous script (for Linux)
net_data_filename = 'net_data.dat';

if ispc == 1                        % If in Windows
    
    evalR('library(igraph);');
    putRdata('P', P);
    putRdata('p', p);
    evalR('G = erdos.renyi.game(P,p);');
    evalR('G = simplify(G,remove.multiple=TRUE,remove.loops=TRUE)');
    evalR('A=get.adjacency(G);');
    Adj = getRdata('A');
    evalR('conn=is.connected(G);');
    conn=getRdata('conn');
    
else                                % If in Linux
    
    % Create python script
    fid = fopen(Erdos_Renyi_filename, 'w');
    
    fprintf(fid, ['import networkx as nx', '\n']);
    fprintf(fid, ['filename = "', net_data_filename, '"' '\n']);
    fprintf(fid, ['P = ', num2str(P) ,'\n']);
    fprintf(fid, ['p = ', num2str(p) ,'\n']);
    fprintf(fid, ['G = nx.generators.erdos_renyi_graph(P,p)', '\n']);    
    fprintf(fid, ['Adj = nx.adj_matrix(G)', '\n']);    
    fprintf(fid, ['if nx.is_connected(G):', '\n']);
    fprintf(fid, ['\tconn = 1', '\n']);
    fprintf(fid, ['else:', '\n']);
    fprintf(fid, ['\tconn = 0', '\n']);
    fprintf(fid, '\n');    
    fprintf(fid, [ 'FILE = open(filename, "w")' ,'\n']);
    fprintf(fid, [ 'FILE.write(''conn = ''+str(int(conn))+''\\n'')','\n']);
    fprintf(fid, [ 'FILE.write(''Adj = \\n'')','\n']);
    fprintf(fid, [ 'for i_1 in range(0, P):','\n']);
    fprintf(fid, [ '\tfor i_2 in range(0, P):','\n']);
    fprintf(fid, [ '\t\tnum = int(Adj[i_1,i_2])','\n']);
    fprintf(fid, [ '\t\tFILE.write(str(num) + '' '')','\n']);
    fprintf(fid, [ '\tFILE.write(''\\n'')','\n']);
    fprintf(fid, '\n');
    fprintf(fid, [ 'FILE.write(''\\n'')','\n']);
    fprintf(fid, [ 'FILE.close()','\n']);
    
    fclose(fid);
    
    system(['python ', Erdos_Renyi_filename]);  % Execute Python script    
    system(['rm ', Erdos_Renyi_filename]);  % Delete script file
    
    % Read data from datafile
    [conn, Adj] = read_Adjacency_and_connectivity(net_data_filename, P);    
    system(['rm ', net_data_filename]); % Delete datafile
    
end

end




function [conn, Adj] = Watts_Strogatz_interface(P, nei, p)

% Python filename (for Linux)
WS_filename = 'Watts_Strogatz.py';

% Filename of output file in case we execute the previous script (for Linux)
net_data_filename = 'net_data.dat';

if ispc == 1                        % If in Windows
    
    evalR('library(igraph);');
    putRdata('P', P);
    putRdata('p', p);
    putRdata('nei', nei);
    evalR('G = watts.strogatz.game(1, P, nei, p);');
    evalR('G = simplify(G,remove.multiple=TRUE,remove.loops=TRUE)');
    evalR('A=get.adjacency(G);');
    Adj = getRdata('A');
    evalR('conn=is.connected(G);');
    conn=getRdata('conn');
    
else                                % If in Linux
    
    % Create python script
    fid = fopen(WS_filename, 'w');
    
    fprintf(fid, ['import networkx as nx', '\n']);
    fprintf(fid, ['filename = "', net_data_filename, '"' '\n']);
    fprintf(fid, ['P = ', num2str(P) ,'\n']);
    fprintf(fid, ['nei = ', num2str(nei) ,'\n']);
    fprintf(fid, ['p = ', num2str(p) ,'\n']);
    
    fprintf(fid, ['G = nx.generators.watts_strogatz_graph(P,nei,p)', '\n']);    
    fprintf(fid, ['Adj = nx.adj_matrix(G)', '\n']);    
    fprintf(fid, ['if nx.is_connected(G):', '\n']);
    fprintf(fid, ['\tconn = 1', '\n']);
    fprintf(fid, ['else:', '\n']);
    fprintf(fid, ['\tconn = 0', '\n']);
    fprintf(fid, '\n');    
    fprintf(fid, [ 'FILE = open(filename, "w")' ,'\n']);
    fprintf(fid, [ 'FILE.write(''conn = ''+str(int(conn))+''\\n'')','\n']);
    fprintf(fid, [ 'FILE.write(''Adj = \\n'')','\n']);
    fprintf(fid, [ 'for i_1 in range(0, P):','\n']);
    fprintf(fid, [ '\tfor i_2 in range(0, P):','\n']);
    fprintf(fid, [ '\t\tnum = int(Adj[i_1,i_2])','\n']);
    fprintf(fid, [ '\t\tFILE.write(str(num) + '' '')','\n']);
    fprintf(fid, [ '\tFILE.write(''\\n'')','\n']);
    fprintf(fid, '\n');
    fprintf(fid, [ 'FILE.write(''\\n'')','\n']);
    fprintf(fid, [ 'FILE.close()','\n']);
    
    fclose(fid);
    
    system(['python ', WS_filename]);  % Execute Python script    
    system(['rm ', WS_filename]);  % Delete script file
    
    % Read data from datafile
    [conn, Adj] = read_Adjacency_and_connectivity(net_data_filename, P);    
    system(['rm ', net_data_filename]); % Delete datafile
    
end

end



function [conn, Adj] = Barabasi_interface(P, m)

% Python filename (for Linux)
Barabasi_filename = 'Barabasi.py';

% Filename of output file in case we execute the previous script (for Linux)
net_data_filename = 'net_data.dat';

if ispc == 1                        % If in Windows
    
    evalR('library(igraph);');
    putRdata('P', P);
    putRdata('m', m);
    evalR('G = barabasi.game(P, power = 1, m, directed=FALSE);');
    evalR('G = simplify(G,remove.multiple=TRUE,remove.loops=TRUE)');
    evalR('A=get.adjacency(G);');
    Adj = getRdata('A');
    evalR('conn=is.connected(G);');
    conn=getRdata('conn');
    
else                                % If in Linux
    
    % Create python script
    fid = fopen(Barabasi_filename, 'w');
    
    fprintf(fid, ['import networkx as nx', '\n']);
    fprintf(fid, ['filename = "', net_data_filename, '"' '\n']);
    fprintf(fid, ['P = ', num2str(P) ,'\n']);
    fprintf(fid, ['m = ', num2str(m) ,'\n']);
    
    fprintf(fid, ['G = nx.generators.barabasi_albert_graph(P,m)', '\n']);    
    fprintf(fid, ['Adj = nx.adj_matrix(G)', '\n']);    
    fprintf(fid, ['if nx.is_connected(G):', '\n']);
    fprintf(fid, ['\tconn = 1', '\n']);
    fprintf(fid, ['else:', '\n']);
    fprintf(fid, ['\tconn = 0', '\n']);
    fprintf(fid, '\n');    
    fprintf(fid, [ 'FILE = open(filename, "w")' ,'\n']);
    fprintf(fid, [ 'FILE.write(''conn = ''+str(int(conn))+''\\n'')','\n']);
    fprintf(fid, [ 'FILE.write(''Adj = \\n'')','\n']);
    fprintf(fid, [ 'for i_1 in range(0, P):','\n']);
    fprintf(fid, [ '\tfor i_2 in range(0, P):','\n']);
    fprintf(fid, [ '\t\tnum = int(Adj[i_1,i_2])','\n']);
    fprintf(fid, [ '\t\tFILE.write(str(num) + '' '')','\n']);
    fprintf(fid, [ '\tFILE.write(''\\n'')','\n']);
    fprintf(fid, '\n');
    fprintf(fid, [ 'FILE.write(''\\n'')','\n']);
    fprintf(fid, [ 'FILE.close()','\n']);
    
    fclose(fid);
    
    system(['python ', Barabasi_filename]);  % Execute Python script    
    system(['rm ', Barabasi_filename]);  % Delete script file
    
    % Read data from datafile
    [conn, Adj] = read_Adjacency_and_connectivity(net_data_filename, P);    
    system(['rm ', net_data_filename]); % Delete datafile
    
end

end



function [conn, Adj] = Geometric_interface(P, d)

% Python filename (for Linux)
Geo_filename = 'Geometric.py';

% Filename of output file in case we execute the previous script (for Linux)
net_data_filename = 'net_data.dat';

if ispc == 1                        % If in Windows
    
    evalR('library(igraph);');
    putRdata('P', P);
    putRdata('d', d);
    evalR('G = grg.game(P, d);');
    evalR('G = simplify(G,remove.multiple=TRUE,remove.loops=TRUE)');
    evalR('A=get.adjacency(G);');
    Adj = getRdata('A');
    evalR('conn=is.connected(G);');
    conn=getRdata('conn');
    
else                                % If in Linux
    
    % Create python script
    fid = fopen(Geo_filename, 'w');
    
    fprintf(fid, ['import networkx as nx', '\n']);
    fprintf(fid, ['filename = "', net_data_filename, '"' '\n']);
    fprintf(fid, ['P = ', num2str(P) ,'\n']);
    fprintf(fid, ['d = ', num2str(d) ,'\n']);
    
    fprintf(fid, ['G = nx.random_geometric_graph(P,d)', '\n']);    
    fprintf(fid, ['Adj = nx.adj_matrix(G)', '\n']);    
    fprintf(fid, ['if nx.is_connected(G):', '\n']);
    fprintf(fid, ['\tconn = 1', '\n']);
    fprintf(fid, ['else:', '\n']);
    fprintf(fid, ['\tconn = 0', '\n']);
    fprintf(fid, '\n');    
    fprintf(fid, [ 'FILE = open(filename, "w")' ,'\n']);
    fprintf(fid, [ 'FILE.write(''conn = ''+str(int(conn))+''\\n'')','\n']);
    fprintf(fid, [ 'FILE.write(''Adj = \\n'')','\n']);
    fprintf(fid, [ 'for i_1 in range(0, P):','\n']);
    fprintf(fid, [ '\tfor i_2 in range(0, P):','\n']);
    fprintf(fid, [ '\t\tnum = int(Adj[i_1,i_2])','\n']);
    fprintf(fid, [ '\t\tFILE.write(str(num) + '' '')','\n']);
    fprintf(fid, [ '\tFILE.write(''\\n'')','\n']);
    fprintf(fid, '\n');
    fprintf(fid, [ 'FILE.write(''\\n'')','\n']);
    fprintf(fid, [ 'FILE.close()','\n']);
    
    fclose(fid);
    
    system(['python ', Geo_filename]);  % Execute Python script    
    system(['rm ', Geo_filename]);  % Delete script file
    
    % Read data from datafile
    [conn, Adj] = read_Adjacency_and_connectivity(net_data_filename, P);    
    system(['rm ', net_data_filename]); % Delete datafile
    
end

end


function [conn, Adj] = Lattice_interface(n1, n2)

% Python filename (for Linux)
% Lattice_filename = 'Lattice.py';

% Filename of output file in case we execute the previous script (for Linux)
% net_data_filename = 'net_data.dat';

P = n1*n2;

n1_new = max(n1,n2);
n2_new = min(n1,n2);

n1 = n1_new;
n2 = n2_new;

d = zeros(1,n1);
d(2) = 1;
D = toeplitz(d);

Adj = zeros(P,P);

Ind = eye(n1);

for i = 1 : n2-1
    indices = 1 + (i-1)*n1 : i*n1;    
    indices2 = 1 + i*n1 : (i+1)*n1;
    
    Adj(indices,indices) = D;
    Adj(indices, indices2) = Ind;
    Adj(indices2, indices) = Ind;
end

i = n2;
indices = 1 + (i-1)*n1 : i*n1;
Adj(indices,indices) = D;
conn = 1;

% if ispc == 1                        % If in Windows
%     
%     evalR('library(igraph);');
%     putRdata('n1', n1);
%     putRdata('n2', n2);
%     evalR('dimvec=rep(0,2);');
%     evalR('dimvec[1] = n1');
%     evalR('dimvec[2] = n2');    
%     evalR('G = graph.lattice(dimvector=dimvec,directed=FALSE);');
%     evalR('A=get.adjacency(G);');
%     Adj = getRdata('A');
%     evalR('conn=is.connected(G);');
%     conn=getRdata('conn');
%     
% else                                % If in Linux
%     
%     % Create python script
%     fid = fopen(Lattice_filename, 'w');
%     
%     fprintf(fid, ['import networkx as nx', '\n']);
%     fprintf(fid, ['filename = "', net_data_filename, '"' '\n']);
%     fprintf(fid, ['n1 = ', num2str(n1) ,'\n']);
%     fprintf(fid, ['n2 = ', num2str(n2) ,'\n']);
%     fprintf(fid, ['P = ', num2str(P) ,'\n']);
%     
%     fprintf(fid, ['G = nx.generators.grid_2d_graph(n1,n2)', '\n']);    
%     fprintf(fid, ['Adj = nx.adj_matrix(G)', '\n']);    
%     fprintf(fid, ['if nx.is_connected(G):', '\n']);
%     fprintf(fid, ['\tconn = 1', '\n']);
%     fprintf(fid, ['else:', '\n']);
%     fprintf(fid, ['\tconn = 0', '\n']);
%     fprintf(fid, '\n');    
%     fprintf(fid, [ 'FILE = open(filename, "w")' ,'\n']);
%     fprintf(fid, [ 'FILE.write(''conn = ''+str(int(conn))+''\\n'')','\n']);
%     fprintf(fid, [ 'FILE.write(''Adj = \\n'')','\n']);
%     fprintf(fid, [ 'for i_1 in range(0, P):','\n']);
%     fprintf(fid, [ '\tfor i_2 in range(0, P):','\n']);
%     fprintf(fid, [ '\t\tnum = int(Adj[i_1,i_2])','\n']);
%     fprintf(fid, [ '\t\tFILE.write(str(num) + '' '')','\n']);
%     fprintf(fid, [ '\tFILE.write(''\\n'')','\n']);
%     fprintf(fid, '\n');
%     fprintf(fid, [ 'FILE.write(''\\n'')','\n']);
%     fprintf(fid, [ 'FILE.close()','\n']);        
%     
%     fclose(fid);
%     
%     system(['python ', Lattice_filename]);  % Execute Python script    
%     system(['rm ', Lattice_filename]);  % Delete script file
%     
%     % Read data from datafile
%     [conn, Adj] = read_Adjacency_and_connectivity(net_data_filename, P);    
%     system(['rm ', net_data_filename]); % Delete datafile
%    
%end

end



function [E] = GN_convertAdj2ListOfNodes(Adj, NN)

% [E] = GN_nconvertAdj2ListOfNodes(Adj, NN)
%
% Converts an adjencency matrix 'Adj' of a graph to a list of nodes 'E' of
% dimensions 2 x NE, where 'NE' is the number of edges. 'NN' is the number
% of nodes.

E = [];

e = 1;

NAdj = Adj + Adj';

for i = 1 : NN
    for j = i+1 : NN
        
        if NAdj(i,j) > 0
            v = [i;j];
            E = [E , v]; %#ok<AGROW>
            e = e + 1;
        end
    end
end

end
