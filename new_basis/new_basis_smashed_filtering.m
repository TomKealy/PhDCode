clear all;
close all;

M = 300;
K = 300:-10:10;
runs = 100;
mse_orth = zeros(size(K, 2), runs);
mse_direct = zeros(size(K, 2), runs);
mse_b = zeros(size(K, 2), runs);
miss_class_runs = zeros(size(K, 2), runs);
miss_class_admm_runs = zeros(size(K, 2), runs);

for i=1:size(K,2)
    edges = [50, 120, 170, 192, 220, 224, 256, 300] ;
    levels = [0,  0 , 200, 0, 0, 0, 0, 0];
    idxs = zeros(1, M)  ;
    idxs(edges(1: end-1)+1) = 1 ;
    g = levels(cumsum(idxs)+1);
    sigma_squared = [1, 10, 100, 1000];
    noise = normrnd(0, 1, [1, M]);
    max_not_in_positions = 0;

    tic
    for j=1:runs
        %     positions = randi(M, [1, 5]);
        %     weights = normrnd(100, sigma_squared(j), [1, 5]);
        %
        %     g = zeros(1,M);
        %     g(positions) = weights;
        
        ge = g + noise;
        
        F = LehmerMatrix(M);
        [L, U] = lu(F);
        I = eye(M);
        D = inv(L);
        
        h = cumsum(g)';
        
        he = cumsum(ge)';
               
        A = normrnd(0, 1/(K(i)), [K(i), M]);
        
        y = A*g';
       
        orth_templates = zeros(M, K(i));
        direct_templates = zeros(M, K(i));
        b = zeros(1, M);
        c = zeros(1, M);
        
        rho = nthroot(max(abs(eig(A'*A))), 3);
        lambda = sqrt(2*log(M));
        alpha = 1;
        
        ops = struct('max_iter', {1000}, ...
        'quiet', {1}, ...
        'abs_tol', {1e-4}, ...
        'rel_tol', {1e-2}....
        );
        
        [z, history] = lasso_admm_1(A*L, y, lambda, 0.5, alpha, ops);
        
        for k=1:M
            orth_templates(k, :) = (A*A')\A*L(k,:)';
            direct_templates(k, :) = A*L(k, :)';
            b(k) = (M/K(i))*(y'*orth_templates(k,:)');
            c(k) = K(i)*(y'*direct_templates(k,:)');
        end
        
        class = zeros(1, M);
        miss_class = zeros(1,M);
        miss_class_admm = zeros(1, M);
        
        ghat = L\b';
        zhat = L*z;
        thresh = 0.5*max(ghat);
        admm_thresh = 0.5*max(zhat);
        for ii=1:M
            if ghat(ii) <= thresh
                ghat(ii) = 0;
            else
                ghat(ii) = 1;
            end
            
            if zhat(ii) <= thresh
                zhat(ii) = 0;
            else 
                zhat(ii) = 1;
            end
            
            if g(ii) == 0
                class(ii) = 0;
            else
                class(ii) = 1;
            end
            
            if ghat(ii) == class(ii)
                miss_class(ii) = 0;
            else
                miss_class(ii) = 1;
            end
            
            if zhat(ii) == class(ii)
                miss_class_admm(ii) = 0;
            else
                miss_class_admm(ii) = 1;
            end
        end
          
        miss_class_runs(i, j) = nnz(miss_class)/M;
        miss_class_admm_runs(i, j) = nnz(miss_class_admm)/M;
        mse_orth(i, j) = norm(b' - h)/norm(h);
        mse_direct(i, j) = norm(c' - h)/norm(h);
        mse_b(i, j) = norm(L\b' - g')/norm(g');
    end
    toc

end

av_b = mean(mse_b, 2);
av_h = mean(mse_orth, 2);
av_dir = mean(mse_direct, 2);
k_over_m = K./M;
plot(av_dir, k_over_m, 'g',  av_b, k_over_m, 'b')
legend('estimate', 'transformed estimate')

figure
plot(mean(miss_class_runs, 2), k_over_m, 'b', mean(miss_class_admm_runs, 2), k_over_m, 'r');


ahat = L'\ghat;



% figure
% plot(ghat)
% 
% figure
% plot(ahat)
% 
% [sort_vals, sort_inds] = sort(abs(ahat), 'descend');
% 
% best_k = sort_inds(1:2);
% 
% for ii=1:M
%     if any(ii == best_k)
%         ahat(ii) = ahat(ii);
%     else
%         ahat(ii) = 0;
%     end
% end
% 
% ghat2 = L'*ahat;
% 
% figure
% plot(ghat2)
