%Script to generate a setof Tx and Rx within the sqaure [0,1]^2, and
%calculate the channel response (Rayleigh fading).

clear all;

num_Tx = 1;
num_Rx = 50;

SNR = 0.1;

Tx_positions = [0.5, 0.5];%rand(num_Tx,2);
Rx_positions = rand(num_Rx,2);

%Generate edges between pairs of within norm of 0.25.
r=0.25
adj_matrix = zeros(num_Rx, num_Rx);
adj_list = cell([1,num_Rx]);
for i=1:num_Rx
    for j=i+1:num_Rx
        dist = norm(Rx_positions(i,:) - Rx_positions(j,:));
        if dist < r;
            adj_matrix(i,j) = 1;
        end
    end
end

for i=1:num_Rx
    I = find( and( adj_matrix(i,:)>0,  adj_matrix(i,:)~=Inf) );
    adj_list{i} = I;
end

figure
gplot(adj_matrix, Rx_positions)

adj_sparse = sparse(adj_matrix);
distances = (graphallshortestpaths(adj_sparse));
dfin = isfinite(distances)

for i=1:num_Rx
    for j=1:num_Rx
        if dfin(i,j) == 1
            ;
        else
            distances(i,j) = 0;
        end
    end
end

% if you don't care which nodes produce the
% longest shortest path, omit long_ind & long_sub
[diameter, long_ind] = max(distances(:));
long_sub = ind2sub(size(distances), long_ind);

%color partition
partition = cell(1,50);
for i=1:num_Rx
    partition{1,i} = i
end

%construct the network struct
Networks = cell(1);

field1 = 'P';
value1 = num_Rx;

field2 = 'Adj';
value2 = adj_matrix;

field3 = 'Type';
value3 = 'Geometric';

field4 = 'parameters'
value4 = 0.5 

field5 = 'Colours';
value5 = num_Rx;

field6 = 'Diameter';
value6 = 10%diameter; % From Diameter and Broadcast Time of Random Geometric Graphs in Arbitrary Dimensions Tobias Friedrich · Thomas Sauerwald · Alexandre Stauffer 

field7 = 'Neighbours';
value7 = adj_list;

field8 = 'Partition';
value8 = partition;

network = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5, field6, value6, field7, value7, field8, value8);

Networks{1,1} = network;

%Tx PSD
W=600; %bandwidth of spectrum
L=201; %number of sub-bands
L0 = (L-1)/2; %re-labelling via shift
B=W/L; %bandwidth per subband
K = 0.02*L; % sparsity

positions = randi(L,[1,K]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1;
figure
plot(Tx_psd)

%deterministic path loss model exp(|x_t - x_r|^alpha)
alpha = 2;
channel_loss = zeros(num_Rx, 1);
dist = zeros(num_Rx, 2);
time_of_flight = zeros(num_Rx,1);
Rx_psd = zeros(num_Rx,L);
S = randsrc(num_Rx, L);
%Loop to calculate channel loss from Tx to ith Rx, and perform mixing
for i=1:num_Rx
    c = 3e8;
    d = Tx_positions(1,:) - Rx_positions(i,:);
    channel_loss(i) = exp(-norm(d,alpha)); 
    dist(i,:) = d;
    Rx_psd(i,:) = channel_loss(i)*Tx_psd + sqrt(SNR)*randn(1,L);
    Rx_psd(i,:) = S(i,:).*Rx_psd(i,:);
end

figure
plot(20*log10(abs(Rx_psd)))

theta = exp(-1i*2*pi/L);
F = theta.^([0:L-1]'*[-L0:L0]);
np = 1:L0;
nn = (-L0):1:-1;
% This is for digital input only. Note that when R -> infinity,
% D then coincides with that of the paper
dn = [   (1-theta.^nn)./(1-theta.^(nn))/(L)      1/L    (1-theta.^np)./(1-theta.^(np))/(L)];
D = diag(dn);
A = S;%*F*D;%diag(channel_loss);
A= conj(A);
A_BP = A;

b = A*Tx_psd';

% % Frame construction
% Q = Rx_psd* Rx_psd';
% % decompose Q to find frame V
% NumDomEigVals= FindNonZeroValues(eig(Q),5e-8);
% [V,d] = eig_r(Q,min(NumDomEigVals,2*L));
% v = V*diag(sqrt(d));
% % N iterations at most, since we force symmetry in the support...
% [u, RecSupp] = RunOMP_Unnormalized(v, A, L, 0, 0.01, true);
% RecSuppSorted = sort(unique(RecSupp));
% if (is_contained(positions,RecSuppSorted)  && (rank(A(:,RecSuppSorted)) == length(RecSuppSorted)))
%     Success= 1;
%     fprintf('Successful support recovery\n');
% else
%     Success = 0;
%     fprintf('Failed support recovery\n');
% end
% 
% A_S = A(:,RecSuppSorted);
% hat_zn = pinv(A_S)*Rx_psd; 
% 
% hat_zt = zeros(size(hat_zn,1),length(Tx_psd));
% for ii = 1:size(hat_zt,1)                     % interpolate (by sinc)
%     hat_zt(ii,:) = interpft(hat_zn(ii,:),L);
% end
% 
% x_rec = zeros(1,length(Tx_psd));
% for ii = 1:size(hat_zt,1)                      % modulate each band to their corresponding carriers
%     x_rec = x_rec+hat_zt(ii,:).*exp(j*2*pi*(RecSuppSorted(ii)));
% end

save('/home/tk12098/Documents/MATLAB/model.mat', 'Networks')

save('/home/tk12098/Documents/MATLAB/psd_data.mat', 'b', 'A_BP', 'Tx_psd')