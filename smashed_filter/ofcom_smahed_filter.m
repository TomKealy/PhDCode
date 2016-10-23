clear all;

data = csvread('tvws_data1.csv');

l = size(data);

data = data(10001:20000, 2);

M = 10000;
K = M/25;

g = data; %reshape(data, [10,1000]);

F = LehmerMatrix(M);
[L, U] = lu(F);
I = eye(M);
D = inv(L);

h = cumsum(g)';

A = normrnd(0, 1/(K), [K, M]);

y = A*g;

direct_templates = zeros(M, K);
c = zeros(10, M);

for ll=1:10
    for k=1:M
        direct_templates(k, :) = A*L(k, :)';
        c(ll, k) = K*(y(:,ll)'*direct_templates(k,:)');
    end
end
class = zeros(1, M);
miss_class = zeros(1, M);

ghat = L\c';
ahat = F\c';
thresh = 0.2*max(ghat);

%plot(1:M, ghat, 'b', 1:M, g, 'r')

% Find k biggest values of ahat

ahat = reshape(ahat, [1,10000]);

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
inds = [1, ahat_big, 10000];
inds = sort(inds);
piece_ttests = zeros(1, num_change_points);

for l=1:num_change_points
    piece = ghat(inds(l):inds(l+1))';
    piece_ttests(l) = ttest(piece, 0, 'Tail', 'right');
end

estimate = zeros(1, 10*M);
threshest = zeros(1, 10*M);

for ii=1:num_change_points
    for jj = inds(ii):inds(ii+1)
        if piece_ttests(ii) == 1
            estimate(jj) = 30; %piece_mean(ii);
        else
            estimate(jj) = -6;
        end
    end
end

errors = zeros(10,M);

g = reshape(g, [1,10000]);

figure
plot(1:10*M, g, 'r', 1:10*M, estimate, 'b')