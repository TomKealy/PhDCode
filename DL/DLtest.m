clear all;

randn('seed', 0);
rand('seed',0);

m=300;

edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [50,  0 , 400, 0, 0, 200, 0, 0];
idxs = zeros(1, m)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1) ;

edges = [50, 120, 170, 192, 220, 234, 256, 300] ;
levels = [100,  0 , 400, 0, 0, 0, 0, 0];
idxs = zeros(1, m)  ;
idxs(edges(1: end-1)+1) = 1 ;
gstar = levels(cumsum(idxs)+1) ;

sigma = 20*10^(-100000\20);
eta = randn(1,m);
noise_sum = sum(eta);

num_examples = 50;
Y = zeros(m, num_examples);

% for ii=1:num_examples
%     decay = randn(1,1);
%     bandwidth = logsig(randn(1,1));
%     rate = randi(50, 1)*10;
%     duration = 10*(rate/11)*1e-3;
%     num_pulses = length(0 : 1/rate : duration)-1;
%     t = 0 : 1/50e3 : 10e-3;
%     d = [0 : 1/rate : duration ; decay.^(0:num_pulses)]';
%     Y(:,ii) = pulstran(t, d,'gauspuls', 15e3, bandwidth);
% end

%csY = Y(:, 1:20);

S = randn(100, 300);
Template = zeros(300, 50);
for ii=1:num_examples
   rr = logsig(randn(1,1));
   if rr < 0.2
       Y(:,ii)  = g;
   else
       Y(:,ii) = gstar;
   end
   Template(:, ii) = gstar;
end

b = S*Y;% + 100*randn(300,num_examples);

% data = csvread('tvws_data1.csv');
% 
% l = size(data);
% start = 100001;
% F = LehmerMatrix(M);
% [L, U] = lu(F);
% I = eye(M);
% D = inv(L);
% 
% chunklength = 10000;
% Y = zeros(10000, 50);
% 
% for ii=1:50
%     
%     Y(:,ii) = data(start:start+(chunklength-1), 2);
%     figure
%     plot(chunk)
%     start = start + chunklength + 3000;
%     estimate = smashed_filt_estimate(chunk, L, F);
%     h = figure;
%     plot(1:M, chunk, 'b', 1:M, estimate, 'r')
%     string = strcat('OFCOM_Chunking', int2str(ii));
%     saveas(h, string, 'epsc')
% end
% 
% S = randn(4000, 10000);
% b = S*Y;
F = LehmerMatrix(300);
[L, U] = lu(F);
[x, D, A] = SparseDL_admm(Y, b, S, 10, 1, Template', L);
%[xbar, Dbar] = CompressiveDL_admm(b, S, D, Y, 1);

% d = A*Y;
% 
% lambda = sqrt(2*log(501));
% 
% max_eig = max(abs(eig((A*D)'*(A*D))));
% 
% rho = nthroot(1/max_eig, 3);
% [xhat, history] = lasso(A*D, d, lambda, rho, 1.0, x, 1000);
% norm(D*xhat-Y)/norm(Y)
% 
% max_eig1 = max(abs(eig((S*D)'*(S*D))));
% 
% rho1 = nthroot(1/max_eig, 3);
% [xhat1, history1] = lasso(S*D, d, lambda, rho1, 1.0, x, 1000);
% norm(D*xhat1-Y)/norm(Y)
% 
% l = length(history.objval)
% 
% figure
% semilogy(1:l, history.objval, 'b', 1:l, history1.objval, 'r')
% legend('optimal', 'random');
% 
% figure
% plot(A'*A)
% title('Optimal Sensing Matrix')
% 
% figure
% plot(S'*S)
% title('Random Sensing Matrix')