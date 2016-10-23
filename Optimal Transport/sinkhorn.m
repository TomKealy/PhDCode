%Parameters
gamma = 10; %reg for entropy
maxiter = 100; % maxiter
map = colormap(gray);
 
 
N = 100; % size
x = linspace(0,1,N)';%spatial coordinate
 
% marginals
p = exp(-(x-0.2).^2*10^2) + exp(-abs(x-0.4)*20);p=p./sum(p); %for colums sum
q = exp(-(x-0.8).^2*10^2);q = q./sum(q); % for row sum
 
[i,j] = meshgrid(1:N);
M = exp(-(i-j).^2/gamma); % exp(-cost/gamma)
 
% intialize u and v
u = ones(N,1);v = ones(N,1);
 
% Sinkhorn-Knopp
% iteratively scale rows and columns
for k = 1:maxiter
    % update u and v
    u = p./(M*v);
    v = q./(M'*u);
    % assemble pi (only for illustration purposes)
    pi = diag(v)*M*diag(u);
    % display pi (with marginals on top and to the left)
    imagesc([p'/max(p) 0;pi/max(pi(:)) q/max(q)])
    colormap(1-map)
    drawnow
end