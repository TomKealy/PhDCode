clear all;

M = 10000;
K = M/2;

edges =  [0, 500, 1500, 2500, 3500, 4500, 5000, 5200, 6000, 7800. 8400, 9000, 10000];
levels = [0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1);

%g = data;

F = LehmerMatrix(M);
[L, U] = lu(F);
I = eye(M);
D = inv(L);

h = cumsum(g)';
runs = 3;
scores = zeros(1, M);

pvals = [0.01, 0.05, 0.1, 0.2, 0.25, 0.3,  0.5, 0.75, 0.8, 0.9, 0.99];
tprs = zeros(runs, length(pvals));
fprs = zeros(runs, length(pvals));

for run = 1:runs
for rr=1:length(pvals)
    A = normrnd(0, 1/(K), [K, M]);
    
    y = A*g';
    
    direct_templates = zeros(M, K);
    c = zeros(1, M);
    
    for k=1:M
        direct_templates(k, :) = A*L(k, :)';
        c(k) = K*(y'*direct_templates(k,:)');
    end
    
    class = zeros(1, M);
    miss_class = zeros(1, M);
    
    ghat = L\c';
    ahat = F\c';
    thresh = 0.2*max(ghat);
    
    % Find k biggest values of ahat
    
    [sorted_ahat, inds] = sort(abs(ahat), 'descend');
    
    num_change_points = 80;
    
    ahat_big = inds(1:num_change_points);
    
    for ii = 1:M
        if ismember(ii, ahat_big);
            ii;
        else
            ahat(ii) = 0;
        end
    end
    
    % take the average of ghat between these vaules
    piece_mean = zeros(1, num_change_points);
    piece_var = zeros(1,num_change_points);
    threshths = zeros(1,num_change_points);
    inds = [1, ahat_big', 300];
    inds = sort(inds);
    piece_ttests = zeros(1, num_change_points);
    ps = zeros(1, num_change_points);
    
    for l=1:num_change_points
        piece = ghat(inds(l):inds(l+1))';
        piece_mean(l) = mean(piece);
        piece_var(l) = var(piece,1);
        threshths(l) = (1.96*piece_var(l))/sqrt(length(piece));
        [piece_ttests(l), ps(l)] = ttest(piece, 0, 'Tail', 'right', 'Alpha', pvals(rr));
    end
    
    estimate = zeros(1, M);
    threshest = zeros(1, M);
    
    for ii=1:num_change_points
        for jj = inds(ii):inds(ii+1)
            if piece_ttests(ii) == 1
                estimate(jj) = 1; %piece_mean(ii);
            else
                estimate(jj) = 0;
            end
        end
    end

    
    tpr = 0;
    fpr = 0;
    
    for kk=1:M
        if estimate(kk) == 1 && g(kk) == 1
            tpr = tpr + 1;
        end
        if estimate(kk) == 1 && g(kk) == 0
            fpr = fpr + 1;
        end
        if estimate(kk) == g(kk)
            scores(kk) = scores(kk) + 1;
        end
    end
    num_ones = nnz(g);
    num_zero = M-num_ones;
    tprs(run, rr) = tpr/num_ones;
    fprs(run, rr) = fpr/num_zero;
end
end
[X, Y] = perfcurve(g, scores, 1);
plot(X, Y);
csvwrite('scores_run1.dat', [tprs;fprs]);