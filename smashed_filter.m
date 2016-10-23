clear all;
close all;

N = 300;
M = 200;

% edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
% levels = [200,  0 , 400, 0, 0, 200, 0, 0];
% idxs = zeros(1, M)  ;
% idxs(edges(1: end-1)+1) = 1 ;
% g = levels(cumsum(idxs)+1);
 
sigma_squared = [1, 10, 100, 1000];
noise = normrnd(0, 1, [1, N]);
runs = 1;% size(sigma_squared,2);

mse = zeros(1,runs);

tic
for j=1:runs
    positions = randi(N, [1, 3]);
    weights = normrnd(100, sigma_squared(1), [1, 1]);
    
    g = zeros(1,N);
    g(positions) = weights;
    
    I = eye(N);

    A = normrnd(0, 1/M, [M, N]);
    
    %time domain measurements
    y = A*g';
       
    templates = zeros(M, N);
    b = zeros(1,N);
    
    for i=1:N
        templates(:, i) = (A*A')\A*I(i,:)';
        b(i) = (y'*templates(:,i));
    end
        
    mse(j) = norm(b-g)/norm(g);
    
    [sortVals, sortInds] = sort(b, 'descend');
    c = sortInds(1:1);
    
    Ind = max(abs(b));
      
    figure
    plot(1:N, abs(b), 'ro', 1:N, g', 'b')
    str = sprintf('Compressive Filtering');
    title(str)
    legend('estimate', 'original')
    
    [U, S, V] = svd(A);
    
    EigDecomp = V'*(S)'*(S)*V;
    
%     %orthogonalised
%     
%     orth_templates = zeros(M, N);
%     b = zeros(1,N);
%     for i=1:N
%         orth_templates(:, i) = (A*A')\A*I(i,:)';
%         c(i) = (y'*orth_templates(:,i));
%     end
%     
%     figure
%     plot(1:N, abs(c), 'ro', 1:N, g', 'b')
%     str = sprintf('Compressive Filtering');
%     title(str)
%         
end
toc