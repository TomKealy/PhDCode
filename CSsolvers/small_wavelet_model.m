clear all;
close all;

M = 256;
K = 200;

edges = [50, 120, 170, 192, 220, 224, 256] ;
levels = [0,  0 , 400, 0, 0, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
Tx_psd = levels(cumsum(idxs)+1)' ;
npsd = Tx_psd;
sigma = 1;
noise = normrnd(0, 1, [1,M]);
Tx_psd = Tx_psd + noise';

% positions = randi(M,[1,5])
% Tx_psd = zeros(1,M)'; %Tx PSD
% Tx_psd(positions) = 1;

%Tx_psd = H(2,:)' - 0.5*H(3,:)' + 0.5*H(4,:)' - H(5,:)' + H(6,:)' - 0.5*H(7,:)' + H(8,:)' - 0.25*H(9,:)' - 1.1*H(3,:)' + 0.5*H(51,:)' - 0.25*H(128,:)' - 1.2*H(35,:)' + H(3,:)' - 0.1*H(201,:)';

% data = csvread('tvws_data1.csv');
% 
% data_length = 2^13;
% x = data(1:data_length,2);
% Tx_psd = x;
% 
% M = length(Tx_psd);
% K = 0.75*M;

H = haarmtx(M); % ConstructHaarWaveletTransformationMatrix(M) ;

% % %d i f f e r e n t i a t i o n matrix
n(1:M-1) = -1;
D = diag(n, -1) + speye(M) ;

% %K Random measurements%
S = rand(K, M);
A_BP = S*H;
%s = spline(1:M, Tx_psd, 1:0.5:M);
%w = wavedec(Tx_psd,8,'haar');       
b = A_BP*Tx_psd;
eg = max(abs(eig((A_BP)'*(A_BP))));
rho = nthroot(1/eg,3);
lambda = sqrt(2*log(M));
relaxation_parameter = 0.7;

[z0_r, history] = lasso_admm_1(S, b, lambda, rho, relaxation_parameter);

initsigma2 = std(b)^2/1e2;
[weights,used,sigma2,errbars,basis,selected,alpha,lambdas] = FastLaplace(S, b, initsigma2, 1e-8, []);

y0 = z0_r;

figure
plot(1:M, H*Tx_psd, 'b', 1:M, y0, 'm')

x_Lap = zeros(M,1);  x_Lap(used) = weights;
err = zeros(M,1); err(used) = errbars;
s_Lap = length(used);
E_Lap = norm(H*npsd-x_Lap)/norm(npsd);
legend('true', 'estimate')

non_zeros = find(H*npsd);

figure(1), subplot(2,1,1), stem(H*npsd), title('Original Signal'), axis([1 M -2 2]);
subplot(2,1,2),stem(x_Lap), axis([1 M -2 2]);, title(sprintf('Reconstruction from %d measurements, error = %g',M,E_Lap));

%partition = [4, 1, 2, 4, 3, 9, 1, 2, 21, 7, 53, 1, 148];

%[x history] = group_lasso(A_BP, b, lambda, partition, rho, 1.0);
