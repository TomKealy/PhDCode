%Nodes = [2 4 8 16 32 64 128 256 512 1024];
Nodes = [1e1 5e1 1e2 2e2 5e2 7e2 1e3 2e3];
parameters = cell(1,1);
type = 'Geometric';
Filename_saving = 'Nets_Geometric_10_2000_TEST.mat';

Num_Netws = length(Nodes);
Networks = cell(Num_Netws,1);

for i_num_net = 1 : Num_Netws
    P = Nodes(i_num_net);
    parameters{1} = sqrt(log(P)/P);
    Networks{i_num_net} = GenerateNetwork(Nodes(i_num_net), type, parameters);
    Networks{i_num_net}
end
save(Filename_saving, 'Networks');














