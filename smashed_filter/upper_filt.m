clear all;
close all;

M = 300;
K = 250; %300:-10:10;
runs = 100;
mse_orth = zeros(size(K, 2), runs);
mse_direct = zeros(size(K, 2), runs);
mse_b = zeros(size(K, 2), runs);
miss_class_runs = zeros(size(K, 2), runs);

edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
levels = [0,  0 , 400, 0, 0, 0, 500, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1);
sigma_squared = [1, 10, 100, 1000];
noise = normrnd(0, 1, [1, M]);
max_not_in_positions = 0;
g = g;% + noise;
tic

%     positions = randi(M, [1, 5]);
%     weights = normrnd(100, sigma_squared(j), [1, 5]);
%
%     g = zeros(1,M);
%     g(positions) = weights;

ge = g + noise;

F = LehmerMatrix(M);
[L U] = lu(F);
I = eye(M);
D = inv(L);

h = cumsum(g)';

he = cumsum(ge)';

A = normrnd(0, 1/(K), [K, M]);

y = A*ge';

orth_templates = zeros(M, K);
direct_templates = zeros(M, K);
b = zeros(1, M);
c = zeros(1, M);

for k=1:M
   b(k) =  K*((A'*y)'*L(k,:)');
   c(k) = g*L(k,:)';
end

template = 100*ones(1, 10);

z = A'*y;

thr = 0.5*max(z)

z = z(z>thr);

out_compressed = filter(template, 1, K*A'*y);
out_uncompressed = filter(template, 1, ge);

thresh = 0.99;

% Compute normalizing factor
u = template*template';

n = 1:300;

% Find matches
matches = n(out_compressed>thresh*u);

% Plot the results
plot(n , out_compressed ,'b', n(matches), out_compressed(matches), 'ro');

% Print the results to the console
display(matches);