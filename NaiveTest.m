clear all
close all 

M = 300;
K = 300;

% % %d i f f e r e n t i a t i o n matrix
n(1:M-1) = -1;
D = diag(n , -1) + speye(M) ;

edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [50,  0 , 400, 0, 0, 200, 0, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1) ;
sigma_squared = 10;
noise = sigma_squared*randn(1, M);

ge = g + noise;

F = LehmerMatrix(M);
[L U] = lu(F);

h = cumsum(g)';

he = cumsum(ge)';

%h = Fa;

A = rand(200,M);

snr = norm(g)/norm(noise);

y = A*he;

eg = max(abs(eig((A'*A))));
rho = nthroot(1/eg,3);
lambdas = 1:ceil(sqrt(log(M))):ceil(sqrt(sigma_squared*2*log(M)));
%lambda = 5*;
relaxation_parameter = 1.0;
path = zeros(size(lambdas,2),M);

for i =1:size(lambdas,2);
    [ghat, history] = lasso_admm_1(A*L, y, 50*sqrt(2*sigma_squared*log(M)), 0.5, relaxation_parameter);
    path(i,:) = ghat;
end
 
figure
plot(1:M, ghat, 'm', 1:M, g, 'b')
str = sprintf('New basis, sigma^2 = %f', sigma_squared);
title(str)

% ghat = L'*ahat;
% 
% figure
% plot(1:M, ghat, 'm', 1:M, g, 'b')
% str = sprintf('New basis, sigma^2 = %f', sigma_squared);
% title(str)